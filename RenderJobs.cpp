#include "gui/jobs/RenderJobs.h"
#include "gravity/BarnesHut.h"
#include "gravity/Moments.h"
#include "gui/Factory.h"
#include "gui/Project.h"
#include "gui/objects/Camera.h"
#include "gui/objects/Movie.h"
#include "gui/objects/PaletteEntry.h"
#include "quantities/Attractor.h"
#include "run/IRun.h"
#include "run/VirtualSettings.h"
#include "run/jobs/IoJobs.h"
#include "system/Factory.h"
#include "system/Timer.h"

#ifdef SPH_USE_VDB
#include <openvdb/openvdb.h>
#endif


NAMESPACE_SPH_BEGIN

//-----------------------------------------------------------------------------------------------------------
// AnimationJob
//-----------------------------------------------------------------------------------------------------------

static RegisterEnum<AnimationType> sAnimation({
    //动画类型
    { AnimationType::SINGLE_FRAME, "single_frame", "Renders only single frame." }, //单帧渲染
    { AnimationType::FILE_SEQUENCE,
        "file_sequence",
        "Make animation from saved files." }, //文件序列：从保存的文件中制作动画
});

static RegisterEnum<RenderColorizerId> sColorizers({
    //着色器ID
    { RenderColorizerId::VELOCITY, "velocity", "Particle velocities" },               //粒子速度
    { RenderColorizerId::ENERGY, "energy", "Specific internal energy" },              //比内能
    { RenderColorizerId::TEMPERATURE, "temperature", "Temperature" },                 //温度
    { RenderColorizerId::DENSITY, "density", "Density" },                             //密度
    { RenderColorizerId::DAMAGE, "damage", "Damage" },                                //损伤
    { RenderColorizerId::GRAVITY, "gravity", "Gravitational acceleration" },          //重力加速度
    { RenderColorizerId::COMPONENT_ID, "component", "Index of connected component" }, //连接组件索引
    { RenderColorizerId::MATERIAL_ID, "material", "Index of material" },              //材料索引
    { RenderColorizerId::BEAUTY, "beauty", "Beauty" },                                //美观
});

static Palette getPaletteFromProject(Project& project, RenderColorizerId id) { //从项目中获取调色板
    AutoPtr<IColorizer> colorizer; //自动指针,用于获取颜色着色器
    if (id == RenderColorizerId::GRAVITY) {
        colorizer = Factory::getColorizer(project, ColorizerId::ACCELERATION); //获取加速度颜色着色器
    } else {
        colorizer = Factory::getColorizer(project, ColorizerId(id)); //否则,获取相应的着色器ID颜色着色器
    }
    return colorizer->getPalette().value(); //返回着色器的调色板
}

AnimationJob::AnimationJob(const String& name)//构造函数，初始化了任务名称、动画类型、着色器、调色板等信息
    : IImageJob(name) {
    animationType = EnumWrapper(AnimationType::SINGLE_FRAME);//设置默认动画类型为单帧
    colorizerId = EnumWrapper(RenderColorizerId::BEAUTY);//设置默认渲染项为美观

    BeautyColorizer colorizer;
    paletteEntry = ExtraEntry(makeAuto<PaletteEntry>(colorizer.getPalette().value()));

    sequence.units = EnumWrapper(UnitEnum::SI);
}

VirtualSettings AnimationJob::
    getSettings() { //获取虚拟设置,它添加了输出目录和文件掩码的选项,并设置了相应的工具提示
    VirtualSettings connector;               ////虚拟连接器
    addGenericCategory(connector, instName); //添加通用类别

    VirtualSettings::Category& outputCat = connector.addCategory("Output"); //输出类别
    outputCat
        .connect("Directory", "directory", directory)         //连接目录
        .setPathType(IVirtualEntry::PathType::DIRECTORY)      //设置路径类型为目录
        .setTooltip("Directory where the images are saved."); //添加提示：保存图像的目录。
    outputCat
        .connect("File mask", "file_mask", fileMask) //连接文件掩码
        .setTooltip(
            "File mask of the created images. Can contain wildcard %d, which is replaced with the number of "
            "the saved image"); //添加提示：创建镜像的文件掩码。可以包含通配符%d，它被替换为保存图像的编号。

    auto particleEnabler = [this] { //粒子渲染启用器
        return gui.get<RendererEnum>(GuiSettingsId::RENDERER) == RendererEnum::PARTICLE;
    };
    auto raymarcherEnabler = [this] { //光线追踪渲染启用器
        return gui.get<RendererEnum>(GuiSettingsId::RENDERER) == RendererEnum::RAYMARCHER;
    };
    auto surfaceEnabler = [this] { //表面渲染启用器
        const RendererEnum type = gui.get<RendererEnum>(GuiSettingsId::RENDERER);
        return type == RendererEnum::RAYMARCHER || type == RendererEnum::MESH;
    };
    auto volumeEnabler = [this] { //体积渲染启用器
        return gui.get<RendererEnum>(GuiSettingsId::RENDERER) == RendererEnum::VOLUME;
    };
    auto raytraceEnabler = [this] { //光线跟踪渲染启动器：光线追踪或体积
        const RendererEnum type = gui.get<RendererEnum>(GuiSettingsId::RENDERER);
        return type == RendererEnum::RAYMARCHER || type == RendererEnum::VOLUME;
    };

    VirtualSettings::Category& rendererCat =
        connector.addCategory("Rendering"); //添加渲染器相关的设置选项,包括渲染器、物理数量和调色板等
    rendererCat.connect<EnumWrapper>("Renderer", gui, GuiSettingsId::RENDERER); //连接渲染器
    rendererCat
        .connect("Quantity", "quantity", colorizerId) //连接数量
        .setTooltip(
            "Physical quantity used to assign values to particles."); //添加提示：用来给粒子数赋值的物理量
    rendererCat
        .connect("Palette", "palette", paletteEntry) //连接调色板
        .setTooltip(
            "Color palette assigning colors to quantity values.") //添加提示：将颜色分配给数量值的调色板
        .setFallback([this] {                                     //设置回退
            // handle backward compatibility 处理向后兼容性
            RawPtr<PaletteEntry> entry = dynamicCast<PaletteEntry>(paletteEntry.getEntry());
            Palette palette = getPaletteFromProject(Project::getInstance(), RenderColorizerId(colorizerId));
            entry->setPalette(palette.subsample(8));
        });
    rendererCat
        .connect("Include surface gravity", "surface_gravity", addSurfaceGravity) //连接表面引力
        .setEnabler([this] {
            return RenderColorizerId(colorizerId) == RenderColorizerId::GRAVITY;
        }) //设置重力加速度渲染启用器
        .setTooltip("Include the surface gravity of the particle itself."); //添加提示：包括粒子本身的表面引力
    rendererCat
        .connect("Include attractors", "attractor_gravity", addAttractorGravity) //连接引力
        .setEnabler([this] {
            return RenderColorizerId(colorizerId) == RenderColorizerId::GRAVITY;
        })                                                   //设置重力加速度渲染启用器
        .setTooltip("Include the gravity from attractors."); //添加提示：包括外部的引力
    rendererCat.connect<bool>("Transparent background", "transparent", transparentBackground); //连接透明背景
    rendererCat.connect<EnumWrapper>("Color mapping", gui, GuiSettingsId::COLORMAP_TYPE); //连接颜色映射
    rendererCat
        .connect<Float>("Logarithmic factor", gui, GuiSettingsId::COLORMAP_LOGARITHMIC_FACTOR) //连接对数因子
        .setEnabler(
            [&] { return gui.get<ColorMapEnum>(GuiSettingsId::COLORMAP_TYPE) == ColorMapEnum::LOGARITHMIC; });
    rendererCat
        .connect<Float>("Bloom intensity", gui, GuiSettingsId::BLOOM_INTENSITY) //连接泛光强度
        .setEnabler(raytraceEnabler);                                           //设置光线跟踪启动器
    rendererCat
        .connect<Float>("Bloom radius [%]", gui, GuiSettingsId::BLOOM_RADIUS) //连接泛光半径
        .setUnits(0.01f)                                                      //单位
        .setEnabler(raytraceEnabler);
    rendererCat
        .connect<Float>("Particle radius", gui, GuiSettingsId::PARTICLE_RADIUS) //连接粒子半径
        .setEnabler(particleEnabler);                                           //连接粒子启动器
    rendererCat.connect<bool>("Antialiasing", gui, GuiSettingsId::ANTIALIASED)
        .setEnabler(particleEnabler);                                    //连接抗锯齿
    rendererCat.connect<bool>("Show key", gui, GuiSettingsId::SHOW_KEY); //连接显示键
    rendererCat
        .connect<int>("Interation count", gui, GuiSettingsId::RAYTRACE_ITERATION_LIMIT) //连接迭代计数
        .setEnabler(raytraceEnabler);
    rendererCat.connect<Float>("Surface level", gui, GuiSettingsId::SURFACE_LEVEL)
        .setEnabler(surfaceEnabler); //连接表面级别
    rendererCat
        .connect<Vector>("Sun position", gui, GuiSettingsId::SURFACE_SUN_POSITION) //连接太阳位置
        .setEnabler(raytraceEnabler);
    rendererCat
        .connect<Float>("Sunlight intensity", gui, GuiSettingsId::SURFACE_SUN_INTENSITY) //连接阳光强度
        .setEnabler(raytraceEnabler);
    rendererCat
        .connect<Float>("Ambient intensity", gui, GuiSettingsId::SURFACE_AMBIENT) //连接环境光强度
        .setEnabler(raytraceEnabler);
    rendererCat
        .connect<Float>("Surface emission", gui, GuiSettingsId::SURFACE_EMISSION) //连接表面发射
        .setEnabler(raymarcherEnabler);
    rendererCat.connect<EnumWrapper>("BRDF", gui, GuiSettingsId::RAYTRACE_BRDF)
        .setEnabler(raymarcherEnabler); //连接BRDF：双向反射分布函数
    rendererCat
        .connect<Float>("Smoothing factor", gui, GuiSettingsId::RAYTRACE_SMOOTH_FACTOR) //连接平滑因子
        .setEnabler(raymarcherEnabler);
    rendererCat
        .connect<bool>("Render as spheres", gui, GuiSettingsId::RAYTRACE_SPHERES) //连接球体渲染
        .setEnabler(raymarcherEnabler);
    rendererCat
        .connect<bool>("Enable shadows", gui, GuiSettingsId::RAYTRACE_SHADOWS) //连接阴影
        .setEnabler(raymarcherEnabler);
    rendererCat
        .connect<Float>("Medium emission [km^-1]", gui, GuiSettingsId::VOLUME_EMISSION) //连接介质发射
        .setUnits(1.e-3_f)
        .setEnabler(volumeEnabler);
    rendererCat
        .connect<Float>("Medium absorption [km^-1]", gui, GuiSettingsId::VOLUME_ABSORPTION) //连接介质吸收[
        .setUnits(1.e-3_f)
        .setEnabler(volumeEnabler);
    rendererCat
        .connect<Float>("Medium scattering [km^-1]", gui, GuiSettingsId::VOLUME_SCATTERING) //连接介质散射
        .setUnits(1.e-3_f)
        .setEnabler(volumeEnabler);
    rendererCat.connect<Float>(
        "Lensing magnitude", gui, GuiSettingsId::RAYTRACE_LENSING_MAGNITUDE); //连接透镜缩放
    rendererCat
        .connect<bool>("Reduce noise", gui, GuiSettingsId::REDUCE_LOWFREQUENCY_NOISE) //连接减少噪音
        .setEnabler(volumeEnabler);

    VirtualSettings::Category& textureCat = connector.addCategory("Texture paths"); //添加纹理路径类别
    textureCat
        .connect<Path>("Background", gui, GuiSettingsId::RAYTRACE_HDRI) //连接背景
        .setEnabler([this] {                                            //设置启动器：体积或光追
            const RendererEnum id = gui.get<RendererEnum>(GuiSettingsId::RENDERER);
            return id == RendererEnum::VOLUME || id == RendererEnum::RAYMARCHER;
        })
        .setPathType(IVirtualEntry::PathType::INPUT_FILE); //设置路径类型为输入文件

    auto sequenceEnabler = [this] {
        return AnimationType(animationType) == AnimationType::FILE_SEQUENCE;
    }; //序列启用器

    VirtualSettings::Category& animationCat = connector.addCategory("Animation");         //添加动画类别
    animationCat.connect<EnumWrapper>("Animation type", "animation_type", animationType); //连接动画类型
    animationCat
        .connect<Path>("First file", "first_file", sequence.firstFile) //连接序列的第一文件
        .setPathType(IVirtualEntry::PathType::INPUT_FILE)              //设置路径类型为输入文件
        .setFileFormats(getInputFormats())                             //设置文件格式
        .setEnabler(sequenceEnabler);                                  //设置序列启动器
    animationCat.connect<EnumWrapper>("Unit system", "units", sequence.units)
        .setEnabler(sequenceEnabler); //连接序列单位系统
    animationCat
        .connect("Interpolated frames", "extra_frames", extraFrames) //连接插值帧数
        .setEnabler(sequenceEnabler)
        .setTooltip(
            "Sets the number of extra frames added between each two state files."); //添加提示：设置每两个状态文件之间添加的额外帧数。

    return connector;
}

class GravityColorizer : public TypedColorizer<Float> { //重力着色器
private:
    SharedPtr<IScheduler> scheduler; //调度器
    BarnesHut gravity;               //重力计算器
    Array<Float> acc;                //加速度数组
    Float G;                         //引力常数
    bool addSurfaceGravity;          //是否添加表面引力
    bool addAttractorGravity;        //是否添加外部引力

public:
    GravityColorizer(const SharedPtr<IScheduler>& scheduler, //构造函数
        const Palette& palette,
        const Float G,
        const bool addSurfaceGravity,
        const bool addAttractorGravity)
        : TypedColorizer<Float>(QuantityId::POSITION, palette)
        , scheduler(scheduler)
        , gravity(0.8_f, MultipoleOrder::OCTUPOLE, SolidSphereKernel{}, 25, 50, G)
        , G(G)
        , addSurfaceGravity(addSurfaceGravity)
        , addAttractorGravity(addAttractorGravity) {}

    virtual void initialize(const Storage& storage, const RefEnum UNUSED(ref)) override { //初始化重写
        acc.resize(storage.getParticleCnt()); //调整加速度数组大小
        acc.fill(0._f);

        // gravitation acceleration from other particles 来自其他粒子的引力加速度
        gravity.build(*scheduler, storage);

        Array<Vector> dv(storage.getParticleCnt()); //速度差分数组
        dv.fill(Vector(0._f));
        Statistics stats;                               //统计
        gravity.evalSelfGravity(*scheduler, dv, stats); //自引力估计
        if (addAttractorGravity) {                      //如果添加内部引力
            Array<Attractor> attractors = viewToArray(storage.getAttractors());
            gravity.evalAttractors(*scheduler, attractors, dv); //外引力估计
        }
        for (Size i = 0; i < dv.size(); ++i) { //求加速度
            acc[i] = getLength(dv[i]);
        }

        if (addSurfaceGravity) { //如果添加外部引力
            // add surface gravity of each particle
            ArrayView<const Float> m = storage.getValue<Float>(QuantityId::MASS);       //质量视图
            ArrayView<const Vector> r = storage.getValue<Vector>(QuantityId::POSITION); //位置视图
            for (Size i = 0; i < r.size(); ++i) {                                       //求加速度
                acc[i] += G * m[i] / sqr(r[i][H]);
            }
        }
    }

    virtual bool isInitialized() const override { //检查是否初始化
        return !acc.empty();
    }

    virtual Rgba evalColor(const Size idx) const override { //估计颜色
        return palette(acc[idx]);
    }

    virtual Optional<Vector> evalVector(const Size UNUSED(idx)) const override { //估计矢量
        return NOTHING;
    }

    virtual String name() const override { //返回名字
        // needs to 'pretend' to be acceleration to work with palette accessor in IR
        // 需要假装是加速度以与IR中的调色板访问器一起工作
        return "Acceleration";
    }
};

RenderParams AnimationJob::getRenderParams(const GuiSettings& gui) const { //获取渲染参数
    SharedPtr<CameraData> camera = this->getInput<CameraData>("camera");   //相机数据
    RenderParams params;                                                   //渲染参数
    params.camera = camera->camera->clone();                               //设置相机
    params.tracker = std::move(camera->tracker);                           //设置跟踪器
    GuiSettings paramGui = gui;                                            // GUI设置
    paramGui.addEntries(camera->overrides);                                //添加相机覆盖
    params.initialize(paramGui);                                           //添加相机覆盖
    return params;
}

class AnimationRenderOutput : public IRenderOutput { //动画渲染输出
private:
    IRunCallbacks& callbacks; //调用回调
    IRenderer& renderer;      //渲染器
    Size iterationCnt;        //迭代次数

    Timer timer; //定时器
    Size iteration = 0;

public:
    AnimationRenderOutput(IRunCallbacks& callbacks, IRenderer& renderer, const Size iterationCnt) //构造函数
        : callbacks(callbacks)
        , renderer(renderer)
        , iterationCnt(iterationCnt) {}

    virtual void update(const Bitmap<Rgba>& bitmap, Array<Label>&& labels, const bool isFinal) override {
        this->update(bitmap.clone(), std::move(labels), isFinal);
    }

    virtual void update(Bitmap<Rgba>&& bitmap,
        Array<Label>&& labels,
        const bool UNUSED(isFinal)) override { //更新设置用户数据、统计、相对进度和时钟时间,调用时间步回调等
        SharedPtr<AnimationFrame> frame = makeShared<AnimationFrame>(); //动画帧
        frame->bitmap = std::move(bitmap);                              //位图
        frame->labels = std::move(labels);                              //标签
        Storage storage;                                                //存储
        storage.setUserData(frame);                                     //用户数据

        Statistics stats;                                                              //统计
        stats.set(StatisticsId::RELATIVE_PROGRESS, Float(++iteration) / iterationCnt); //设置相对进度
        stats.set(StatisticsId::WALLCLOCK_TIME, int(timer.elapsed(TimerUnit::MILLISECOND))); //设置时钟时间
        callbacks.onTimeStep(storage, stats); //时间步长回调

        if (callbacks.shouldAbortRun()) { //中止运行
            renderer.cancelRender();      //取消渲染
        }
    }
};

void AnimationJob::evaluate(const RunSettings& global, IRunCallbacks& callbacks) { //估计
    if (directory.empty()) {                                                       //如果输出目录为空
        throw InvalidSetup(                                                        //抛出无效设置异常
            L"No output directory specified. Please set the output directory to where you want to save the " //没有指定输出目录。请将输出目录设置为要保存渲染图像的位置。
            "rendered images.");
    }

    /// \todo maybe also work with a copy of Gui ? 也许还可以使用Gui的副本
    gui.set(GuiSettingsId::BACKGROUND_COLOR,
        Rgba(0.f, 0.f, 0.f, transparentBackground ? 0.f : 1.f)); //设置背景颜色
    gui.set(GuiSettingsId::RAYTRACE_SUBSAMPLING, 0);             //设置光线跟踪子采样
    int iterLimit = 1;                                           //迭代限制
    if (gui.get<RendererEnum>(GuiSettingsId::RENDERER) != RendererEnum::PARTICLE) { //如果渲染器不是粒子
        iterLimit = gui.get<int>(GuiSettingsId::RAYTRACE_ITERATION_LIMIT); //迭代限制为光线跟踪迭代限制
    }

    SharedPtr<IScheduler> scheduler = Factory::getScheduler(global);    //设置调度器为全局调度器
    AutoPtr<IRenderer> renderer = Factory::getRenderer(scheduler, gui); //设置渲染器为全局gui渲染器
    RawPtr<IRenderer> rendererPtr = renderer.get();                     //获取渲染器指针

    RenderParams params = this->getRenderParams(gui);           //获取渲染参数
    AutoPtr<IColorizer> colorizer = this->getColorizer(global); //设置着色器为全局着色器

    int firstIndex = 0;
    if (AnimationType(animationType) == AnimationType::FILE_SEQUENCE) { //如果动画类型是文件序列
        Optional<Size> sequenceFirstIndex = OutputFile::getDumpIdx(
            sequence.firstFile);  //获取 sequence.firstFile 的可选输出文件 sequenceFirstIndex
        if (sequenceFirstIndex) { //如果 sequenceFirstIndex 存在,
            firstIndex = sequenceFirstIndex.value(); //将firstIndex的值设置为sequenceFirstIndex的值。
        }
    }
    OutputFile paths(directory / Path(fileMask), firstIndex);            //输出文件路径
    SharedPtr<CameraData> camera = this->getInput<CameraData>("camera"); //获取相机数据
    Movie movie(camera->overrides,
        std::move(renderer),
        std::move(colorizer),
        std::move(params),
        extraFrames,
        paths); //相机覆盖：移动渲染器，移动着色器，移动渲染参数

    switch (AnimationType(animationType)) {
    case AnimationType::SINGLE_FRAME: { //如果动画类型为单帧渲染
        SharedPtr<ParticleData> data = this->getInput<ParticleData>("particles"); //获取粒子数据
        AnimationRenderOutput output(callbacks, *rendererPtr, iterLimit);         //动画渲染输出
        movie.render(std::move(data->storage), std::move(data->stats), output);   //执行渲染
        break;
    }
    case AnimationType::FILE_SEQUENCE: {                                   //如果动画类型为文件序列
        FlatMap<Size, Path> fileMap = getFileSequence(sequence.firstFile); //设置平面映射
        if (fileMap.empty()) {
            throw InvalidSetup("No files to render."); //抛出无效设置：无渲染文件。
        }

        const Size iterationCnt = iterLimit * fileMap.size() * (extraFrames + 1); //迭代计数
        AnimationRenderOutput output(callbacks, *rendererPtr, iterationCnt);      //动画渲染输出
        AutoPtr<IInput> input = Factory::getInput(sequence.firstFile);            //输入文件
        for (auto& element : fileMap) {                                           //对每个文件
            Storage frame;                                                        //帧
            Statistics stats;                                                     //统计
            const Outcome result = input->load(element.value(), frame, stats);    //结果
            if (!result) {                                                        //如果结果不成功
                /// \todo how to report this? (don't do modal dialog) 如何报告?(不要做模态对话框)
            }

            if (callbacks.shouldAbortRun()) { //如果应中止运行
                break;                        //中断
            }

            movie.render(std::move(frame), std::move(stats), output); //执行渲染
        }
        break;
    }
    default:
        NOT_IMPLEMENTED;
    }
}

class RenderPreview : public IRenderPreview { //渲染预览
private:
    RenderParams params;           //渲染参数
    AutoPtr<IRenderer> renderer;   //渲染器
    AutoPtr<IColorizer> colorizer; //着色器
    SharedPtr<ParticleData> data;  //粒子数据
    std::atomic_bool cancelled;    //取消状态

    bool rendererDirty = true;  //渲染器脏
    bool colorizerDirty = true; //着色器脏

public:
    RenderPreview(RenderParams&& params, //构造函数
        AutoPtr<IRenderer>&& renderer,
        AutoPtr<IColorizer>&& colorizer,
        const SharedPtr<ParticleData>& data)
        : params(std::move(params))
        , renderer(std::move(renderer))
        , colorizer(std::move(colorizer))
        , data(data)
        , cancelled(false) {}

    virtual void render(const Pixel resolution, IRenderOutput& output) override { //重写渲染
        cancelled = false;                                                        //未取消

        // lazy init
        if (colorizerDirty) {                                    //延迟初始化
            colorizer->initialize(data->storage, RefEnum::WEAK); //初始化着色器
            colorizerDirty = false;                              //关闭着色器脏
            rendererDirty = true;
        }
        if (cancelled) { //如果取消了，则中断
            return;
        }
        if (rendererDirty) {                                                 //如果渲染器脏
            renderer->initialize(data->storage, *colorizer, *params.camera); //初始化渲染器
            rendererDirty = false;                                           //关闭渲染器
        }
        if (cancelled) {
            return;
        }

        Pixel size = params.camera->getSize();                                //获取尺寸
        size = correctAspectRatio(resolution, float(size.x) / float(size.y)); //纠正宽高比
        params.camera->resize(size);                                          //调整相机尺寸
        Statistics dummy;                                                     //虚拟统计
        renderer->render(params, dummy, output);                              //渲染
    }

    virtual void update(RenderParams&& newParams) override { //重写参数更新
        AutoPtr<ICamera> camera = std::move(params.camera);
        params = std::move(newParams); //参数覆盖
        params.camera = std::move(camera);
    }

    virtual void update(AutoPtr<ICamera>&& newCamera) override { //重写相机更新
        params.camera = std::move(newCamera);                    //相机覆盖
    }

    virtual void update(AutoPtr<IColorizer>&& newColorizer) override { //重写着色器更新
        colorizer = std::move(newColorizer);                           //着色器覆盖
        colorizerDirty = true;                                         //开启着色器脏
    }

    virtual void update(AutoPtr<IRenderer>&& newRenderer) override { //重写渲染器更新
        renderer = std::move(newRenderer);                           //渲染器覆盖
        rendererDirty = true;                                        //开启渲染器脏
    }

    virtual void update(Palette&& palette) override { //重写调色板更新
        colorizer->setPalette(palette);               //着色器获得新调色板
        renderer->setColorizer(*colorizer);           //渲染器获得着色器
    }

    virtual void cancel() override { //重写取消
        cancelled = true;
        renderer->cancelRender();
    }

private:
    Pixel correctAspectRatio(const Pixel resolution, const float aspect) const { //纠正宽高比
        const float current = float(resolution.x) / float(resolution.y);         //获取当前宽高比
        if (current > aspect) { //如果现在的宽高比大于预期的
            return Pixel(resolution.x * aspect / current, resolution.y); //调整宽高比
        } else {                                                         //小于
            return Pixel(resolution.x, resolution.y * current / aspect);
        }
    }
};

AutoPtr<IRenderPreview> AnimationJob::getRenderPreview(const RunSettings& global) const { //获取渲染预览
    if (AnimationType(animationType) != AnimationType::SINGLE_FRAME) { //如果动画类型不是单帧
        throw InvalidSetup("Only enabled for single-frame renders"); //抛出无效设置异常：只适用于单帧渲染。
    }

    if (!inputs.contains("particles")) {               //如果输入不包含粒子
        throw InvalidSetup("Particles not connected"); //抛出无效设置异常：粒子未连接
    }

    RenderParams params = this->getRenderParams(); //获取渲染参数

    AutoPtr<IColorizer> colorizer = this->getColorizer(global); //设置着色器为全局着色器。
    if (!colorizer) {                                           //如果无着色器
        throw InvalidSetup("No quantity selected"); //抛出无效设置异常：未选择数量。
    }
    AutoPtr<IRenderer> renderer = this->getRenderer(global); //设置渲染器为全局渲染器

    SharedPtr<ParticleData> data = this->getInput<ParticleData>("particles"); //获取粒子数据

    return makeAuto<RenderPreview>(
        std::move(params), std::move(renderer), std::move(colorizer), data); //返回渲染预览
}

AutoPtr<IColorizer> AnimationJob::getColorizer(const RunSettings& global) const { //获取着色器
    CHECK_FUNCTION(CheckFunction::NO_THROW);
    RenderColorizerId renderId(colorizerId);
    if (renderId == RenderColorizerId::GRAVITY) {
        SharedPtr<IScheduler> scheduler = Factory::getScheduler(global);
        Float G = Constants::gravity;
        switch (AnimationType(animationType)) {
        case AnimationType::SINGLE_FRAME: {
            SharedPtr<ParticleData> data = this->getInput<ParticleData>("particles");
            if (data->overrides.has(RunSettingsId::GRAVITY_CONSTANT)) {
                G = data->overrides.get<Float>(RunSettingsId::GRAVITY_CONSTANT);
            }
            break;
        }
        case AnimationType::FILE_SEQUENCE:
            G = getGravityConstant(UnitEnum(sequence.units));
            break;
        default:
            NOT_IMPLEMENTED;
        }
        return makeAuto<GravityColorizer>(
            scheduler, this->getPalette(), G, addSurfaceGravity, addAttractorGravity);
    } else {
        AutoPtr<IColorizer> colorizer = Factory::getColorizer(gui, ColorizerId(renderId));
        colorizer->setPalette(this->getPalette());
        return colorizer;
    }
}

AutoPtr<IRenderer> AnimationJob::getRenderer(const RunSettings& global) const { //获取渲染器
    SharedPtr<IScheduler> scheduler = Factory::getScheduler(global);
    GuiSettings previewGui = gui;
    previewGui.set(GuiSettingsId::RAYTRACE_SUBSAMPLING, 4);
    previewGui.set(GuiSettingsId::BACKGROUND_COLOR, Rgba(0.f, 0.f, 0.f, transparentBackground ? 0.f : 1.f));
    AutoPtr<IRenderer> renderer = Factory::getRenderer(scheduler, previewGui);
    return renderer;
}

Palette AnimationJob::getPalette() const { //获取调色板
    RawPtr<PaletteEntry> palette = dynamicCast<PaletteEntry>(paletteEntry.getEntry());
    return palette->getPalette();
}

RenderParams AnimationJob::getRenderParams() const { //获取渲染参数
    GuiSettings previewGui = gui;
    previewGui.set(GuiSettingsId::SHOW_KEY, false);
    previewGui.set(GuiSettingsId::BACKGROUND_COLOR, Rgba(0.f, 0.f, 0.f, transparentBackground ? 0.f : 1.f));
    return this->getRenderParams(previewGui);
}

JobRegistrar sRegisterAnimation( //作业注册表
    "render animation",
    "animation",
    "rendering",
    [](const String& name) { return makeAuto<AnimationJob>(name); },
    "Renders an image or a sequence of images from given particle input(s)");

//-----------------------------------------------------------------------------------------------------------
// VdbJob
//-----------------------------------------------------------------------------------------------------------

#ifdef SPH_USE_VDB


INLINE openvdb::Vec3f vectorToVec3f(const Vector& v) {
    return openvdb::Vec3f(v[X], v[Y], v[Z]);
}

INLINE Vector worldToGrid(const Vector& r, const Box& box, const Indices& dims) {
    return (r - box.lower()) / box.size() * Vector(dims);
}

INLINE Vector gridToWorld(const Vector& r, const Box& box, const Indices& dims) {
    return r * box.size() / Vector(dims) + box.lower();
}

Tuple<Indices, Indices> getParticleBox(const Vector& r, const Box& box, const Indices& dims) {
    const Vector from = worldToGrid(r - Vector(2._f * r[H]), box, dims);
    const Vector to = worldToGrid(r + Vector(2._f * r[H]), box, dims);
    const Indices fromIdxs(ceil(from[X]), ceil(from[Y]), ceil(from[Z]));
    const Indices toIdxs(floor(to[X]), floor(to[Y]), floor(to[Z]));
    return { max(fromIdxs, Indices(0._f)), min(toIdxs, dims - Indices(1)) };
}

VirtualSettings VdbJob::getSettings() {
    VirtualSettings connector;
    addGenericCategory(connector, instName);

    VirtualSettings::Category& gridCat = connector.addCategory("Grid parameters");
    gridCat.connect("Grid start [km]", "grid_start", gridStart)
        .setUnits(1.e3_f)
        .setTooltip("Sets the lower bound of the bounding box.");
    gridCat.connect("Grid end [km]", "grid_end", gridEnd)
        .setUnits(1.e3_f)
        .setTooltip("Sets the upper bound of the bounding box.");
    gridCat.connect("Resolution power", "power", dimPower)
        .setTooltip("Defines resolution of the grid. The number of voxels in one dimension is 2^power.");
    gridCat.connect("Surface level", "surface_level", surfaceLevel).setTooltip("Iso-value of the surface.");

    VirtualSettings::Category& inputCat = connector.addCategory("File sequence");
    inputCat.connect("Enable", "enable_sequence", sequence.enabled);
    inputCat.connect("First file", "first_file", sequence.firstFile)
        .setPathType(IVirtualEntry::PathType::INPUT_FILE)
        .setFileFormats(getInputFormats())
        .setEnabler([this] { return sequence.enabled; });

    VirtualSettings::Category& outputCat = connector.addCategory("Output");
    outputCat.connect("VDB File", "file", path)
        .setPathType(IVirtualEntry::PathType::OUTPUT_FILE)
        .setFileFormats({ { "OpenVDB grid file", "vdb" } })
        .setEnabler([this] { return !sequence.enabled; });

    return connector;
}

void VdbJob::evaluate(const RunSettings& global, IRunCallbacks& callbacks) {
    openvdb::initialize();
    auto deinit = finally([] { openvdb::uninitialize(); });

    if (sequence.enabled) {
        FlatMap<Size, Path> fileMap = getFileSequence(sequence.firstFile);
        if (fileMap.empty()) {
            throw InvalidSetup("No files to render.");
        }
        const Size firstKey = fileMap.begin()->key();

        AutoPtr<IInput> input = Factory::getInput(sequence.firstFile);
        for (auto& element : fileMap) {
            Storage storage;
            Statistics stats;
            const Outcome result = input->load(element.value(), storage, stats);
            if (!result) {
                /// \todo how to report this? (don't do modal dialog)
            }

            Path outputPath = element.value();
            outputPath.replaceExtension("vdb");
            this->generate(storage, global, outputPath);

            /// \todo deduplicate with AnimationJob
            stats.set(StatisticsId::RELATIVE_PROGRESS, Float(element.key() - firstKey) / fileMap.size());
            if (element.key() == firstKey) {
                callbacks.onSetUp(storage, stats);
            }
            callbacks.onTimeStep(storage, stats);

            if (callbacks.shouldAbortRun()) {
                break;
            }
        }
    } else {
        Storage& storage = getInput<ParticleData>("particles")->storage;
        this->generate(storage, global, path);
    }
}

void VdbJob::generate(Storage& storage, const RunSettings& global, const Path& outputPath) {
    using namespace openvdb;

    FloatGrid::Ptr colorField = FloatGrid::create(-surfaceLevel);
    Vec3SGrid::Ptr velocityField = Vec3SGrid::create(vectorToVec3f(Vector(0._f)));
    FloatGrid::Ptr energyField = FloatGrid::create(0._f);

    colorField->setName("Density");
    velocityField->setName("Velocity");
    energyField->setName("Emission");

    ArrayView<const Vector> r = storage.getValue<Vector>(QuantityId::POSITION);
    ArrayView<const Vector> v = storage.getDt<Vector>(QuantityId::POSITION);
    ArrayView<const Float> m = storage.getValue<Float>(QuantityId::MASS);
    ArrayView<const Float> u = storage.getValue<Float>(QuantityId::ENERGY);
    ArrayView<const Float> rho = storage.getValue<Float>(QuantityId::DENSITY);

    const Box box(gridStart, gridEnd);
    const Size gridSize = 1 << dimPower;
    const Indices gridIdxs(gridSize);

    LutKernel<3> kernel = Factory::getKernel<3>(global);

    typename FloatGrid::Accessor colorAccessor = colorField->getAccessor();
    typename Vec3SGrid::Accessor velocityAccessor = velocityField->getAccessor();
    typename FloatGrid::Accessor energyAccessor = energyField->getAccessor();
    for (Size i = 0; i < r.size(); ++i) {
        Indices from, to;
        tieToTuple(from, to) = getParticleBox(r[i], box, gridIdxs);
        Float rho_i;
        if (storage.getMaterialCnt() > 0) {
            rho_i = storage.getMaterialOfParticle(i)->getParam<Float>(BodySettingsId::DENSITY);
        } else {
            rho_i = rho[i];
        }
        for (int x = from[X]; x <= to[X]; ++x) {
            for (int y = from[Y]; y <= to[Y]; ++y) {
                for (int z = from[Z]; z <= to[Z]; ++z) {
                    const Indices idxs(x, y, z);
                    const Vector pos = gridToWorld(idxs, box, gridIdxs);
                    const Float w = kernel.value(r[i] - pos, r[i][H]);
                    const Float c = m[i] / rho_i * w;

                    const Coord coord(x, y, z);
                    colorAccessor.modifyValue(coord, [c](float& color) { color += c; });
                    energyAccessor.modifyValue(coord, [&u, c, i](float& energy) { energy += c * u[i]; });
                    velocityAccessor.modifyValue(
                        coord, [&v, c, i](Vec3f& velocity) { velocity += c * vectorToVec3f(v[i]); });
                }
            }
        }
    }

    for (FloatGrid::ValueOnIter iter = colorField->beginValueOn(); iter; ++iter) {
        const Coord coord = iter.getCoord();
        const float c = *iter;
        if (c > 0) {
            energyAccessor.modifyValue(coord, [c](float& energy) { energy /= c; });
            velocityAccessor.modifyValue(coord, [c](Vec3f& velocity) { velocity /= c; });
        }
        iter.setValue(c - surfaceLevel);
    }

    GridPtrVec vdbGrids;
    vdbGrids.push_back(colorField);
    vdbGrids.push_back(velocityField);
    vdbGrids.push_back(energyField);

    Path vdbPath = outputPath;
    vdbPath.replaceExtension("vdb");
    io::File vdbFile(vdbPath.string().toAscii().cstr());
    vdbFile.write(vdbGrids);
    vdbFile.close();
}

JobRegistrar sRegisterVdb(
    "save VDB grid",
    "grid",
    "rendering",
    [](const String& name) { return makeAuto<VdbJob>(name); },
    "Converts the particle data into a volumetric grid in OpenVDB format.");

#endif

NAMESPACE_SPH_END
