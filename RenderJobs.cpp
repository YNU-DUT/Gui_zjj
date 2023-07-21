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
    //��������
    { AnimationType::SINGLE_FRAME, "single_frame", "Renders only single frame." }, //��֡��Ⱦ
    { AnimationType::FILE_SEQUENCE,
        "file_sequence",
        "Make animation from saved files." }, //�ļ����У��ӱ�����ļ�����������
});

static RegisterEnum<RenderColorizerId> sColorizers({
    //��ɫ��ID
    { RenderColorizerId::VELOCITY, "velocity", "Particle velocities" },               //�����ٶ�
    { RenderColorizerId::ENERGY, "energy", "Specific internal energy" },              //������
    { RenderColorizerId::TEMPERATURE, "temperature", "Temperature" },                 //�¶�
    { RenderColorizerId::DENSITY, "density", "Density" },                             //�ܶ�
    { RenderColorizerId::DAMAGE, "damage", "Damage" },                                //����
    { RenderColorizerId::GRAVITY, "gravity", "Gravitational acceleration" },          //�������ٶ�
    { RenderColorizerId::COMPONENT_ID, "component", "Index of connected component" }, //�����������
    { RenderColorizerId::MATERIAL_ID, "material", "Index of material" },              //��������
    { RenderColorizerId::BEAUTY, "beauty", "Beauty" },                                //����
});

static Palette getPaletteFromProject(Project& project, RenderColorizerId id) { //����Ŀ�л�ȡ��ɫ��
    AutoPtr<IColorizer> colorizer; //�Զ�ָ��,���ڻ�ȡ��ɫ��ɫ��
    if (id == RenderColorizerId::GRAVITY) {
        colorizer = Factory::getColorizer(project, ColorizerId::ACCELERATION); //��ȡ���ٶ���ɫ��ɫ��
    } else {
        colorizer = Factory::getColorizer(project, ColorizerId(id)); //����,��ȡ��Ӧ����ɫ��ID��ɫ��ɫ��
    }
    return colorizer->getPalette().value(); //������ɫ���ĵ�ɫ��
}

AnimationJob::AnimationJob(const String& name)//���캯������ʼ�����������ơ��������͡���ɫ������ɫ�����Ϣ
    : IImageJob(name) {
    animationType = EnumWrapper(AnimationType::SINGLE_FRAME);//����Ĭ�϶�������Ϊ��֡
    colorizerId = EnumWrapper(RenderColorizerId::BEAUTY);//����Ĭ����Ⱦ��Ϊ����

    BeautyColorizer colorizer;
    paletteEntry = ExtraEntry(makeAuto<PaletteEntry>(colorizer.getPalette().value()));

    sequence.units = EnumWrapper(UnitEnum::SI);
}

VirtualSettings AnimationJob::
    getSettings() { //��ȡ��������,����������Ŀ¼���ļ������ѡ��,����������Ӧ�Ĺ�����ʾ
    VirtualSettings connector;               ////����������
    addGenericCategory(connector, instName); //���ͨ�����

    VirtualSettings::Category& outputCat = connector.addCategory("Output"); //������
    outputCat
        .connect("Directory", "directory", directory)         //����Ŀ¼
        .setPathType(IVirtualEntry::PathType::DIRECTORY)      //����·������ΪĿ¼
        .setTooltip("Directory where the images are saved."); //�����ʾ������ͼ���Ŀ¼��
    outputCat
        .connect("File mask", "file_mask", fileMask) //�����ļ�����
        .setTooltip(
            "File mask of the created images. Can contain wildcard %d, which is replaced with the number of "
            "the saved image"); //�����ʾ������������ļ����롣���԰���ͨ���%d�������滻Ϊ����ͼ��ı�š�

    auto particleEnabler = [this] { //������Ⱦ������
        return gui.get<RendererEnum>(GuiSettingsId::RENDERER) == RendererEnum::PARTICLE;
    };
    auto raymarcherEnabler = [this] { //����׷����Ⱦ������
        return gui.get<RendererEnum>(GuiSettingsId::RENDERER) == RendererEnum::RAYMARCHER;
    };
    auto surfaceEnabler = [this] { //������Ⱦ������
        const RendererEnum type = gui.get<RendererEnum>(GuiSettingsId::RENDERER);
        return type == RendererEnum::RAYMARCHER || type == RendererEnum::MESH;
    };
    auto volumeEnabler = [this] { //�����Ⱦ������
        return gui.get<RendererEnum>(GuiSettingsId::RENDERER) == RendererEnum::VOLUME;
    };
    auto raytraceEnabler = [this] { //���߸�����Ⱦ������������׷�ٻ����
        const RendererEnum type = gui.get<RendererEnum>(GuiSettingsId::RENDERER);
        return type == RendererEnum::RAYMARCHER || type == RendererEnum::VOLUME;
    };

    VirtualSettings::Category& rendererCat =
        connector.addCategory("Rendering"); //�����Ⱦ����ص�����ѡ��,������Ⱦ�������������͵�ɫ���
    rendererCat.connect<EnumWrapper>("Renderer", gui, GuiSettingsId::RENDERER); //������Ⱦ��
    rendererCat
        .connect("Quantity", "quantity", colorizerId) //��������
        .setTooltip(
            "Physical quantity used to assign values to particles."); //�����ʾ����������������ֵ��������
    rendererCat
        .connect("Palette", "palette", paletteEntry) //���ӵ�ɫ��
        .setTooltip(
            "Color palette assigning colors to quantity values.") //�����ʾ������ɫ���������ֵ�ĵ�ɫ��
        .setFallback([this] {                                     //���û���
            // handle backward compatibility ������������
            RawPtr<PaletteEntry> entry = dynamicCast<PaletteEntry>(paletteEntry.getEntry());
            Palette palette = getPaletteFromProject(Project::getInstance(), RenderColorizerId(colorizerId));
            entry->setPalette(palette.subsample(8));
        });
    rendererCat
        .connect("Include surface gravity", "surface_gravity", addSurfaceGravity) //���ӱ�������
        .setEnabler([this] {
            return RenderColorizerId(colorizerId) == RenderColorizerId::GRAVITY;
        }) //�����������ٶ���Ⱦ������
        .setTooltip("Include the surface gravity of the particle itself."); //�����ʾ���������ӱ���ı�������
    rendererCat
        .connect("Include attractors", "attractor_gravity", addAttractorGravity) //��������
        .setEnabler([this] {
            return RenderColorizerId(colorizerId) == RenderColorizerId::GRAVITY;
        })                                                   //�����������ٶ���Ⱦ������
        .setTooltip("Include the gravity from attractors."); //�����ʾ�������ⲿ������
    rendererCat.connect<bool>("Transparent background", "transparent", transparentBackground); //����͸������
    rendererCat.connect<EnumWrapper>("Color mapping", gui, GuiSettingsId::COLORMAP_TYPE); //������ɫӳ��
    rendererCat
        .connect<Float>("Logarithmic factor", gui, GuiSettingsId::COLORMAP_LOGARITHMIC_FACTOR) //���Ӷ�������
        .setEnabler(
            [&] { return gui.get<ColorMapEnum>(GuiSettingsId::COLORMAP_TYPE) == ColorMapEnum::LOGARITHMIC; });
    rendererCat
        .connect<Float>("Bloom intensity", gui, GuiSettingsId::BLOOM_INTENSITY) //���ӷ���ǿ��
        .setEnabler(raytraceEnabler);                                           //���ù��߸���������
    rendererCat
        .connect<Float>("Bloom radius [%]", gui, GuiSettingsId::BLOOM_RADIUS) //���ӷ���뾶
        .setUnits(0.01f)                                                      //��λ
        .setEnabler(raytraceEnabler);
    rendererCat
        .connect<Float>("Particle radius", gui, GuiSettingsId::PARTICLE_RADIUS) //�������Ӱ뾶
        .setEnabler(particleEnabler);                                           //��������������
    rendererCat.connect<bool>("Antialiasing", gui, GuiSettingsId::ANTIALIASED)
        .setEnabler(particleEnabler);                                    //���ӿ����
    rendererCat.connect<bool>("Show key", gui, GuiSettingsId::SHOW_KEY); //������ʾ��
    rendererCat
        .connect<int>("Interation count", gui, GuiSettingsId::RAYTRACE_ITERATION_LIMIT) //���ӵ�������
        .setEnabler(raytraceEnabler);
    rendererCat.connect<Float>("Surface level", gui, GuiSettingsId::SURFACE_LEVEL)
        .setEnabler(surfaceEnabler); //���ӱ��漶��
    rendererCat
        .connect<Vector>("Sun position", gui, GuiSettingsId::SURFACE_SUN_POSITION) //����̫��λ��
        .setEnabler(raytraceEnabler);
    rendererCat
        .connect<Float>("Sunlight intensity", gui, GuiSettingsId::SURFACE_SUN_INTENSITY) //��������ǿ��
        .setEnabler(raytraceEnabler);
    rendererCat
        .connect<Float>("Ambient intensity", gui, GuiSettingsId::SURFACE_AMBIENT) //���ӻ�����ǿ��
        .setEnabler(raytraceEnabler);
    rendererCat
        .connect<Float>("Surface emission", gui, GuiSettingsId::SURFACE_EMISSION) //���ӱ��淢��
        .setEnabler(raymarcherEnabler);
    rendererCat.connect<EnumWrapper>("BRDF", gui, GuiSettingsId::RAYTRACE_BRDF)
        .setEnabler(raymarcherEnabler); //����BRDF��˫����ֲ�����
    rendererCat
        .connect<Float>("Smoothing factor", gui, GuiSettingsId::RAYTRACE_SMOOTH_FACTOR) //����ƽ������
        .setEnabler(raymarcherEnabler);
    rendererCat
        .connect<bool>("Render as spheres", gui, GuiSettingsId::RAYTRACE_SPHERES) //����������Ⱦ
        .setEnabler(raymarcherEnabler);
    rendererCat
        .connect<bool>("Enable shadows", gui, GuiSettingsId::RAYTRACE_SHADOWS) //������Ӱ
        .setEnabler(raymarcherEnabler);
    rendererCat
        .connect<Float>("Medium emission [km^-1]", gui, GuiSettingsId::VOLUME_EMISSION) //���ӽ��ʷ���
        .setUnits(1.e-3_f)
        .setEnabler(volumeEnabler);
    rendererCat
        .connect<Float>("Medium absorption [km^-1]", gui, GuiSettingsId::VOLUME_ABSORPTION) //���ӽ�������[
        .setUnits(1.e-3_f)
        .setEnabler(volumeEnabler);
    rendererCat
        .connect<Float>("Medium scattering [km^-1]", gui, GuiSettingsId::VOLUME_SCATTERING) //���ӽ���ɢ��
        .setUnits(1.e-3_f)
        .setEnabler(volumeEnabler);
    rendererCat.connect<Float>(
        "Lensing magnitude", gui, GuiSettingsId::RAYTRACE_LENSING_MAGNITUDE); //����͸������
    rendererCat
        .connect<bool>("Reduce noise", gui, GuiSettingsId::REDUCE_LOWFREQUENCY_NOISE) //���Ӽ�������
        .setEnabler(volumeEnabler);

    VirtualSettings::Category& textureCat = connector.addCategory("Texture paths"); //�������·�����
    textureCat
        .connect<Path>("Background", gui, GuiSettingsId::RAYTRACE_HDRI) //���ӱ���
        .setEnabler([this] {                                            //������������������׷
            const RendererEnum id = gui.get<RendererEnum>(GuiSettingsId::RENDERER);
            return id == RendererEnum::VOLUME || id == RendererEnum::RAYMARCHER;
        })
        .setPathType(IVirtualEntry::PathType::INPUT_FILE); //����·������Ϊ�����ļ�

    auto sequenceEnabler = [this] {
        return AnimationType(animationType) == AnimationType::FILE_SEQUENCE;
    }; //����������

    VirtualSettings::Category& animationCat = connector.addCategory("Animation");         //��Ӷ������
    animationCat.connect<EnumWrapper>("Animation type", "animation_type", animationType); //���Ӷ�������
    animationCat
        .connect<Path>("First file", "first_file", sequence.firstFile) //�������еĵ�һ�ļ�
        .setPathType(IVirtualEntry::PathType::INPUT_FILE)              //����·������Ϊ�����ļ�
        .setFileFormats(getInputFormats())                             //�����ļ���ʽ
        .setEnabler(sequenceEnabler);                                  //��������������
    animationCat.connect<EnumWrapper>("Unit system", "units", sequence.units)
        .setEnabler(sequenceEnabler); //�������е�λϵͳ
    animationCat
        .connect("Interpolated frames", "extra_frames", extraFrames) //���Ӳ�ֵ֡��
        .setEnabler(sequenceEnabler)
        .setTooltip(
            "Sets the number of extra frames added between each two state files."); //�����ʾ������ÿ����״̬�ļ�֮����ӵĶ���֡����

    return connector;
}

class GravityColorizer : public TypedColorizer<Float> { //������ɫ��
private:
    SharedPtr<IScheduler> scheduler; //������
    BarnesHut gravity;               //����������
    Array<Float> acc;                //���ٶ�����
    Float G;                         //��������
    bool addSurfaceGravity;          //�Ƿ���ӱ�������
    bool addAttractorGravity;        //�Ƿ�����ⲿ����

public:
    GravityColorizer(const SharedPtr<IScheduler>& scheduler, //���캯��
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

    virtual void initialize(const Storage& storage, const RefEnum UNUSED(ref)) override { //��ʼ����д
        acc.resize(storage.getParticleCnt()); //�������ٶ������С
        acc.fill(0._f);

        // gravitation acceleration from other particles �����������ӵ��������ٶ�
        gravity.build(*scheduler, storage);

        Array<Vector> dv(storage.getParticleCnt()); //�ٶȲ������
        dv.fill(Vector(0._f));
        Statistics stats;                               //ͳ��
        gravity.evalSelfGravity(*scheduler, dv, stats); //����������
        if (addAttractorGravity) {                      //�������ڲ�����
            Array<Attractor> attractors = viewToArray(storage.getAttractors());
            gravity.evalAttractors(*scheduler, attractors, dv); //����������
        }
        for (Size i = 0; i < dv.size(); ++i) { //����ٶ�
            acc[i] = getLength(dv[i]);
        }

        if (addSurfaceGravity) { //�������ⲿ����
            // add surface gravity of each particle
            ArrayView<const Float> m = storage.getValue<Float>(QuantityId::MASS);       //������ͼ
            ArrayView<const Vector> r = storage.getValue<Vector>(QuantityId::POSITION); //λ����ͼ
            for (Size i = 0; i < r.size(); ++i) {                                       //����ٶ�
                acc[i] += G * m[i] / sqr(r[i][H]);
            }
        }
    }

    virtual bool isInitialized() const override { //����Ƿ��ʼ��
        return !acc.empty();
    }

    virtual Rgba evalColor(const Size idx) const override { //������ɫ
        return palette(acc[idx]);
    }

    virtual Optional<Vector> evalVector(const Size UNUSED(idx)) const override { //����ʸ��
        return NOTHING;
    }

    virtual String name() const override { //��������
        // needs to 'pretend' to be acceleration to work with palette accessor in IR
        // ��Ҫ��װ�Ǽ��ٶ�����IR�еĵ�ɫ�������һ����
        return "Acceleration";
    }
};

RenderParams AnimationJob::getRenderParams(const GuiSettings& gui) const { //��ȡ��Ⱦ����
    SharedPtr<CameraData> camera = this->getInput<CameraData>("camera");   //�������
    RenderParams params;                                                   //��Ⱦ����
    params.camera = camera->camera->clone();                               //�������
    params.tracker = std::move(camera->tracker);                           //���ø�����
    GuiSettings paramGui = gui;                                            // GUI����
    paramGui.addEntries(camera->overrides);                                //����������
    params.initialize(paramGui);                                           //����������
    return params;
}

class AnimationRenderOutput : public IRenderOutput { //������Ⱦ���
private:
    IRunCallbacks& callbacks; //���ûص�
    IRenderer& renderer;      //��Ⱦ��
    Size iterationCnt;        //��������

    Timer timer; //��ʱ��
    Size iteration = 0;

public:
    AnimationRenderOutput(IRunCallbacks& callbacks, IRenderer& renderer, const Size iterationCnt) //���캯��
        : callbacks(callbacks)
        , renderer(renderer)
        , iterationCnt(iterationCnt) {}

    virtual void update(const Bitmap<Rgba>& bitmap, Array<Label>&& labels, const bool isFinal) override {
        this->update(bitmap.clone(), std::move(labels), isFinal);
    }

    virtual void update(Bitmap<Rgba>&& bitmap,
        Array<Label>&& labels,
        const bool UNUSED(isFinal)) override { //���������û����ݡ�ͳ�ơ���Խ��Ⱥ�ʱ��ʱ��,����ʱ�䲽�ص���
        SharedPtr<AnimationFrame> frame = makeShared<AnimationFrame>(); //����֡
        frame->bitmap = std::move(bitmap);                              //λͼ
        frame->labels = std::move(labels);                              //��ǩ
        Storage storage;                                                //�洢
        storage.setUserData(frame);                                     //�û�����

        Statistics stats;                                                              //ͳ��
        stats.set(StatisticsId::RELATIVE_PROGRESS, Float(++iteration) / iterationCnt); //������Խ���
        stats.set(StatisticsId::WALLCLOCK_TIME, int(timer.elapsed(TimerUnit::MILLISECOND))); //����ʱ��ʱ��
        callbacks.onTimeStep(storage, stats); //ʱ�䲽���ص�

        if (callbacks.shouldAbortRun()) { //��ֹ����
            renderer.cancelRender();      //ȡ����Ⱦ
        }
    }
};

void AnimationJob::evaluate(const RunSettings& global, IRunCallbacks& callbacks) { //����
    if (directory.empty()) {                                                       //������Ŀ¼Ϊ��
        throw InvalidSetup(                                                        //�׳���Ч�����쳣
            L"No output directory specified. Please set the output directory to where you want to save the " //û��ָ�����Ŀ¼���뽫���Ŀ¼����ΪҪ������Ⱦͼ���λ�á�
            "rendered images.");
    }

    /// \todo maybe also work with a copy of Gui ? Ҳ������ʹ��Gui�ĸ���
    gui.set(GuiSettingsId::BACKGROUND_COLOR,
        Rgba(0.f, 0.f, 0.f, transparentBackground ? 0.f : 1.f)); //���ñ�����ɫ
    gui.set(GuiSettingsId::RAYTRACE_SUBSAMPLING, 0);             //���ù��߸����Ӳ���
    int iterLimit = 1;                                           //��������
    if (gui.get<RendererEnum>(GuiSettingsId::RENDERER) != RendererEnum::PARTICLE) { //�����Ⱦ����������
        iterLimit = gui.get<int>(GuiSettingsId::RAYTRACE_ITERATION_LIMIT); //��������Ϊ���߸��ٵ�������
    }

    SharedPtr<IScheduler> scheduler = Factory::getScheduler(global);    //���õ�����Ϊȫ�ֵ�����
    AutoPtr<IRenderer> renderer = Factory::getRenderer(scheduler, gui); //������Ⱦ��Ϊȫ��gui��Ⱦ��
    RawPtr<IRenderer> rendererPtr = renderer.get();                     //��ȡ��Ⱦ��ָ��

    RenderParams params = this->getRenderParams(gui);           //��ȡ��Ⱦ����
    AutoPtr<IColorizer> colorizer = this->getColorizer(global); //������ɫ��Ϊȫ����ɫ��

    int firstIndex = 0;
    if (AnimationType(animationType) == AnimationType::FILE_SEQUENCE) { //��������������ļ�����
        Optional<Size> sequenceFirstIndex = OutputFile::getDumpIdx(
            sequence.firstFile);  //��ȡ sequence.firstFile �Ŀ�ѡ����ļ� sequenceFirstIndex
        if (sequenceFirstIndex) { //��� sequenceFirstIndex ����,
            firstIndex = sequenceFirstIndex.value(); //��firstIndex��ֵ����ΪsequenceFirstIndex��ֵ��
        }
    }
    OutputFile paths(directory / Path(fileMask), firstIndex);            //����ļ�·��
    SharedPtr<CameraData> camera = this->getInput<CameraData>("camera"); //��ȡ�������
    Movie movie(camera->overrides,
        std::move(renderer),
        std::move(colorizer),
        std::move(params),
        extraFrames,
        paths); //������ǣ��ƶ���Ⱦ�����ƶ���ɫ�����ƶ���Ⱦ����

    switch (AnimationType(animationType)) {
    case AnimationType::SINGLE_FRAME: { //�����������Ϊ��֡��Ⱦ
        SharedPtr<ParticleData> data = this->getInput<ParticleData>("particles"); //��ȡ��������
        AnimationRenderOutput output(callbacks, *rendererPtr, iterLimit);         //������Ⱦ���
        movie.render(std::move(data->storage), std::move(data->stats), output);   //ִ����Ⱦ
        break;
    }
    case AnimationType::FILE_SEQUENCE: {                                   //�����������Ϊ�ļ�����
        FlatMap<Size, Path> fileMap = getFileSequence(sequence.firstFile); //����ƽ��ӳ��
        if (fileMap.empty()) {
            throw InvalidSetup("No files to render."); //�׳���Ч���ã�����Ⱦ�ļ���
        }

        const Size iterationCnt = iterLimit * fileMap.size() * (extraFrames + 1); //��������
        AnimationRenderOutput output(callbacks, *rendererPtr, iterationCnt);      //������Ⱦ���
        AutoPtr<IInput> input = Factory::getInput(sequence.firstFile);            //�����ļ�
        for (auto& element : fileMap) {                                           //��ÿ���ļ�
            Storage frame;                                                        //֡
            Statistics stats;                                                     //ͳ��
            const Outcome result = input->load(element.value(), frame, stats);    //���
            if (!result) {                                                        //���������ɹ�
                /// \todo how to report this? (don't do modal dialog) ��α���?(��Ҫ��ģ̬�Ի���)
            }

            if (callbacks.shouldAbortRun()) { //���Ӧ��ֹ����
                break;                        //�ж�
            }

            movie.render(std::move(frame), std::move(stats), output); //ִ����Ⱦ
        }
        break;
    }
    default:
        NOT_IMPLEMENTED;
    }
}

class RenderPreview : public IRenderPreview { //��ȾԤ��
private:
    RenderParams params;           //��Ⱦ����
    AutoPtr<IRenderer> renderer;   //��Ⱦ��
    AutoPtr<IColorizer> colorizer; //��ɫ��
    SharedPtr<ParticleData> data;  //��������
    std::atomic_bool cancelled;    //ȡ��״̬

    bool rendererDirty = true;  //��Ⱦ����
    bool colorizerDirty = true; //��ɫ����

public:
    RenderPreview(RenderParams&& params, //���캯��
        AutoPtr<IRenderer>&& renderer,
        AutoPtr<IColorizer>&& colorizer,
        const SharedPtr<ParticleData>& data)
        : params(std::move(params))
        , renderer(std::move(renderer))
        , colorizer(std::move(colorizer))
        , data(data)
        , cancelled(false) {}

    virtual void render(const Pixel resolution, IRenderOutput& output) override { //��д��Ⱦ
        cancelled = false;                                                        //δȡ��

        // lazy init
        if (colorizerDirty) {                                    //�ӳٳ�ʼ��
            colorizer->initialize(data->storage, RefEnum::WEAK); //��ʼ����ɫ��
            colorizerDirty = false;                              //�ر���ɫ����
            rendererDirty = true;
        }
        if (cancelled) { //���ȡ���ˣ����ж�
            return;
        }
        if (rendererDirty) {                                                 //�����Ⱦ����
            renderer->initialize(data->storage, *colorizer, *params.camera); //��ʼ����Ⱦ��
            rendererDirty = false;                                           //�ر���Ⱦ��
        }
        if (cancelled) {
            return;
        }

        Pixel size = params.camera->getSize();                                //��ȡ�ߴ�
        size = correctAspectRatio(resolution, float(size.x) / float(size.y)); //������߱�
        params.camera->resize(size);                                          //��������ߴ�
        Statistics dummy;                                                     //����ͳ��
        renderer->render(params, dummy, output);                              //��Ⱦ
    }

    virtual void update(RenderParams&& newParams) override { //��д��������
        AutoPtr<ICamera> camera = std::move(params.camera);
        params = std::move(newParams); //��������
        params.camera = std::move(camera);
    }

    virtual void update(AutoPtr<ICamera>&& newCamera) override { //��д�������
        params.camera = std::move(newCamera);                    //�������
    }

    virtual void update(AutoPtr<IColorizer>&& newColorizer) override { //��д��ɫ������
        colorizer = std::move(newColorizer);                           //��ɫ������
        colorizerDirty = true;                                         //������ɫ����
    }

    virtual void update(AutoPtr<IRenderer>&& newRenderer) override { //��д��Ⱦ������
        renderer = std::move(newRenderer);                           //��Ⱦ������
        rendererDirty = true;                                        //������Ⱦ����
    }

    virtual void update(Palette&& palette) override { //��д��ɫ�����
        colorizer->setPalette(palette);               //��ɫ������µ�ɫ��
        renderer->setColorizer(*colorizer);           //��Ⱦ�������ɫ��
    }

    virtual void cancel() override { //��дȡ��
        cancelled = true;
        renderer->cancelRender();
    }

private:
    Pixel correctAspectRatio(const Pixel resolution, const float aspect) const { //������߱�
        const float current = float(resolution.x) / float(resolution.y);         //��ȡ��ǰ��߱�
        if (current > aspect) { //������ڵĿ�߱ȴ���Ԥ�ڵ�
            return Pixel(resolution.x * aspect / current, resolution.y); //������߱�
        } else {                                                         //С��
            return Pixel(resolution.x, resolution.y * current / aspect);
        }
    }
};

AutoPtr<IRenderPreview> AnimationJob::getRenderPreview(const RunSettings& global) const { //��ȡ��ȾԤ��
    if (AnimationType(animationType) != AnimationType::SINGLE_FRAME) { //����������Ͳ��ǵ�֡
        throw InvalidSetup("Only enabled for single-frame renders"); //�׳���Ч�����쳣��ֻ�����ڵ�֡��Ⱦ��
    }

    if (!inputs.contains("particles")) {               //������벻��������
        throw InvalidSetup("Particles not connected"); //�׳���Ч�����쳣������δ����
    }

    RenderParams params = this->getRenderParams(); //��ȡ��Ⱦ����

    AutoPtr<IColorizer> colorizer = this->getColorizer(global); //������ɫ��Ϊȫ����ɫ����
    if (!colorizer) {                                           //�������ɫ��
        throw InvalidSetup("No quantity selected"); //�׳���Ч�����쳣��δѡ��������
    }
    AutoPtr<IRenderer> renderer = this->getRenderer(global); //������Ⱦ��Ϊȫ����Ⱦ��

    SharedPtr<ParticleData> data = this->getInput<ParticleData>("particles"); //��ȡ��������

    return makeAuto<RenderPreview>(
        std::move(params), std::move(renderer), std::move(colorizer), data); //������ȾԤ��
}

AutoPtr<IColorizer> AnimationJob::getColorizer(const RunSettings& global) const { //��ȡ��ɫ��
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

AutoPtr<IRenderer> AnimationJob::getRenderer(const RunSettings& global) const { //��ȡ��Ⱦ��
    SharedPtr<IScheduler> scheduler = Factory::getScheduler(global);
    GuiSettings previewGui = gui;
    previewGui.set(GuiSettingsId::RAYTRACE_SUBSAMPLING, 4);
    previewGui.set(GuiSettingsId::BACKGROUND_COLOR, Rgba(0.f, 0.f, 0.f, transparentBackground ? 0.f : 1.f));
    AutoPtr<IRenderer> renderer = Factory::getRenderer(scheduler, previewGui);
    return renderer;
}

Palette AnimationJob::getPalette() const { //��ȡ��ɫ��
    RawPtr<PaletteEntry> palette = dynamicCast<PaletteEntry>(paletteEntry.getEntry());
    return palette->getPalette();
}

RenderParams AnimationJob::getRenderParams() const { //��ȡ��Ⱦ����
    GuiSettings previewGui = gui;
    previewGui.set(GuiSettingsId::SHOW_KEY, false);
    previewGui.set(GuiSettingsId::BACKGROUND_COLOR, Rgba(0.f, 0.f, 0.f, transparentBackground ? 0.f : 1.f));
    return this->getRenderParams(previewGui);
}

JobRegistrar sRegisterAnimation( //��ҵע���
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
