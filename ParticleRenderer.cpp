#include "gui/renderers/ParticleRenderer.h"
#include "gui/Utils.h"
#include "gui/objects/Camera.h"
#include "gui/objects/Color.h"
#include "gui/objects/Colorizer.h"
#include "gui/objects/RenderContext.h"
#include "objects/finders/Order.h"
#include "objects/wrappers/Finally.h"
#include "post/Plot.h"
#include "post/Point.h"
#include "quantities/Attractor.h"
#include "sph/boundary/Boundary.h"
#include "system/Profiler.h"
#include "system/Statistics.h"
#include "thread/CheckFunction.h"

NAMESPACE_SPH_BEGIN
//粒子渲染：用于渲染粒子效果
static void drawVector(IRenderContext& context,//绘制一个带箭头线段表示向量
    const ICamera& camera,
    const Vector& r,
    const Vector& v,
    const Float length) {
    if (getSqrLength(v) == 0._f) {// 如果向量v长度为0,不绘制
        return;
    }
    const Optional<ProjectedPoint> p1 = camera.project(r); // 投影向量起点
    const Optional<ProjectedPoint> p2 = camera.project(r + v);// 投影终点(起点+方向向量)
    if (!p1 || !p2) {//只要p1和p2中有一个无效，则停止绘制。
        return;
    }

    Coords dir = p2->coords - p1->coords;//计算p1指向p2方向的向量dir,
    const Float l = getLength(dir);//获取方向向量长度
    if (l == 0._f) {//如果长度为0，则不绘制
        return;
    }
    dir *= float(length / l);// 将dir缩放到指定长度 length
    const Coords c1 = p1->coords;// 起点坐标
    const Coords c2 = p1->coords + dir;// 终点坐标(起点+方向向量)

    context.setColor(Rgba(1.f, 0.65f, 0.f), ColorFlag::LINE);// 设置线段颜色
    context.setThickness(2.f);// 设置线段粗细
    context.drawLine(c1, c2);// 绘制线段

    // make an arrow 画箭头
    AffineMatrix2 rot = AffineMatrix2::rotate(160._f * DEG_TO_RAD);// 创建一个顺时针旋转160度的变换矩阵
    PlotPoint dp(dir.x, dir.y);// 将方向向量转换为点,以进行矩阵变换
    PlotPoint a1 = rot.transformPoint(dp) * 0.1f;// 对dp点应用旋转矩阵,得到箭头一侧的偏移量
    PlotPoint a2 = rot.transpose().transformPoint(dp) * 0.1f;// 对dp点应用转置旋转矩阵,得到箭头另一侧的偏移量 

    context.drawLine(c2, c2 + Coords(float(a1.x), float(a1.y)));// 在终点出绘制指向a1的线段,形成箭头的一侧
    context.drawLine(c2, c2 + Coords(float(a2.x), float(a2.y)));// 在终点出绘制指向a1的线段,形成箭头的一侧
}

static void drawGrid(IRenderContext& context, const ICamera& camera, const float grid) {//绘制一个与相机对齐的网格背景
    // find (any) direction in the camera plane 在相机平面上找到一个方向
    const Optional<CameraRay> originRay = camera.unproject(Coords(0, 0));// 从屏幕中心发出一条射线
    const Vector dir = getNormalized(originRay->target - originRay->origin);// 计算射线方向向量dir,并归一化
    // 计算与dir垂直的向量perpDir
    Vector perpDir;
    if (dir == Vector(0._f, 0._f, 1._f)) {// 如果dir垂直向上
        perpDir = Vector(1._f, 0._f, 0._f);//则perpDir取(1,0,0)
    } else {
        perpDir = getNormalized(cross(dir, Vector(0._f, 0._f, 1._f)));// 否则取dir与(0,0,1)的叉积作为perpDir
    }

    // find how much is projected grid distance 求出投影的网格距离
    const Coords shifted = camera.project(originRay->origin + grid * perpDir)->coords;//计算原点移动一个网格长度在屏幕上的投影偏移量
    const float dx = getLength(shifted);//取dx的值为shifted，作为垂直网格线单位间距
    const float dy = dx;// dy取与dx相同大小，作为水平网格线单位间距
    const Coords origin = camera.project(Vector(0._f))->coords;// 原点坐标 

    context.setColor(Rgba(0.16f), ColorFlag::LINE);// 设置绘制颜色
    const Pixel size = context.size();// 设置绘制颜色
    for (float x = origin.x; x < size.x; x += dx) {// 在x方向上，从origin.x开始,步长dx,绘制垂直网格线
        context.drawLine(Coords(x, 0), Coords(x, size.y));
    }
    for (float x = origin.x - dx; x >= 0; x -= dx) {//以同样方法在反方向上循环绘制
        context.drawLine(Coords(x, 0), Coords(x, size.y));
    }
    for (float y = origin.y; y < size.y; y += dy) {// 在y方向上，从origin.y开始,步长dy,绘制水平网格线
        context.drawLine(Coords(0, y), Coords(size.x, y));
    }
    for (float y = origin.y - dy; y >= 0; y -= dy) {//以同样方法在反方向上循环绘制
        context.drawLine(Coords(0, y), Coords(size.x, y));
    }
}

ParticleRenderer::ParticleRenderer(const GuiSettings& settings) {//ParticleRenderer的构造函数
    grid = float(settings.get<Float>(GuiSettingsId::VIEW_GRID_SIZE)); // 从GuiSettings对象中,获取VIEW_GRID_SIZE参数指定的网格间距值
  // 该值将用于后续绘制网格背景时控制网格线的间距
    shouldContinue = true;// 初始化shouldContinue标志为true
  // 该标志用于控制渲染循环的继续执行
}

static bool isCutOff(const Vector& r, const Optional<float> cutoff, const Vector direction) {//判断一个向量是否在截断半径外
    return cutoff && abs(dot(direction, r)) > cutoff.value();//如果指定了cutoff,且r在direction方向上的投影长度的绝对值大于cutoff,则返回true,表示在截断半径外。
    //否则返回false,表示在截断半径内。
}

constexpr Size GHOST_INDEX = Size(-1);//虚拟粒子索引
constexpr Size ATTRACTOR_INDEX = Size(-2);//吸引子索引

void ParticleRenderer::initialize(const Storage& storage,//粒子渲染器初始化
    const IColorizer& colorizer,
    const ICamera& camera) {
    MEASURE_SCOPE("ParticleRenderer::initialize");// 记录函数运行时间
    cached.idxs.clear();// 清空粒子索引缓存
    cached.positions.clear(); // 清空粒子位置缓存 
    cached.colors.clear(); // 清空粒子颜色缓存
    cached.vectors.clear();// 清空粒子着色器向量缓存

    //遍历粒子数据,筛选在截断范围内的,提前缓存其屏幕空间信息,为后续渲染做准备。
    const Optional<float> cutoff = camera.getCutoff();// 获取相机的截断半径cutoff
    const Vector direction = camera.getFrame().row(2);//获取相机的视角方向direction
    bool hasVectorData = bool(colorizer.evalVector(0));// 判断着色器是否有向量数据
    ArrayView<const Vector> r = storage.getValue<Vector>(QuantityId::POSITION);//获取所有粒子的位置数据r
    for (Size i = 0; i < r.size(); ++i) {//遍历每个粒子r[i]
        const Optional<ProjectedPoint> p = camera.project(r[i]);// 把粒子位置投影到屏幕
        if (p && !isCutOff(r[i], cutoff, direction)) {// 如果在截断半径内则继续
            cached.idxs.push(i);// 缓存粒子索引
            cached.positions.push(r[i]);//缓存粒子位置

            const Rgba color = colorizer.evalColor(i);//计算粒子颜色
            cached.colors.push(color);//缓存粒子颜色

            if (hasVectorData) {//如果着色器有向量数据
                Optional<Vector> v = colorizer.evalVector(i);//获取单个粒子的着色器向量数据
                SPH_ASSERT(v);//断言向量数据不为空
                cached.vectors.push(v.value());//缓存着色器向量数据
            }
        }
    }

    //处理虚拟粒子
    SharedPtr<IStorageUserData> data = storage.getUserData();// 从存储中获取自定义用户数据 
    if (RawPtr<GhostParticlesData> ghosts = dynamicCast<GhostParticlesData>(data.get())) {// 判断是否是GhostParticlesData类型，从而判断是不是虚拟粒子
        for (Size i = 0; i < ghosts->size(); ++i) {// 遍历每个虚拟粒子
            const Vector pos = ghosts->getGhost(i).position;// 获取虚拟粒子位置信息
            if (!isCutOff(pos, cutoff, direction)) {// 如果在相机截断半径内
                cached.idxs.push(GHOST_INDEX);// 将GHOST_INDEX作为粒子索引缓存
                cached.positions.push(pos);//缓存虚拟粒子的位置信息
                cached.colors.push(Rgba::transparent());//将其颜色设置为透明

                if (hasVectorData) {//如果着色器有向量数据
                    cached.vectors.push(Vector(0._f));//将向量数据设置为0
                }
            }
        }
    }

    //处理吸引子
    for (const Attractor& a : storage.getAttractors()) {//获取存储中的所有吸引子并遍历
        if (!isCutOff(a.position, cutoff, direction)) {// 如果在相机截断半径内
            cached.idxs.push(ATTRACTOR_INDEX);// 将ATTRACTOR_INDEX作为粒子索引缓存
            cached.positions.push(setH(a.position, a.radius));//缓存吸引子的位置信息
            cached.colors.push(Rgba::white());//将其颜色设置为白色

            if (hasVectorData) {//如果着色器有向量数据
                cached.vectors.push(Vector(0._f));//将向量数据设置为0
            }
        }
    }

    // sort in z-order 按Z轴顺序排序，以让所有缓存数据保持同步排序，满足渲染的需求
    Order order(cached.positions.size());//创建一个Order对象,其大小为缓存中位置信息的大小
    order.shuffle([this, &direction](Size i, Size j) {// 使用shuffle方法以二元谓词遍历索引 
        // 比较在视角方向上的投影大小
        const Vector r1 = cached.positions[i];
        const Vector r2 = cached.positions[j];
        return dot(direction, r1) > dot(direction, r2);// 比较两个粒子在视角方向上的投影大小
    });
    /// \todo could be changed to AOS to sort only once 可以改成AOS算法以只排序一次
    // 对位置、索引、颜色应用排序结果
    cached.positions = order.apply(cached.positions);
    cached.idxs = order.apply(cached.idxs);
    cached.colors = order.apply(cached.colors);

    cached.cameraDir = direction;// 对位置、索引、颜色应用排序结果

    if (hasVectorData) {// 如果着色器有向量数据
        cached.vectors = order.apply(cached.vectors);// 对向量数据也应用排序
    } else {
        cached.vectors.clear();//否则清空向量数据
    }
}

bool ParticleRenderer::isInitialized() const {//判断是否初始化
    return !cached.positions.empty();//通过判断cached.positions是否为空,来推断cached中的数据是否已被初始化
}

void ParticleRenderer::setColorizer(const IColorizer& colorizer) {//设置着色器
    for (Size i = 0; i < cached.idxs.size(); ++i) {// 遍历缓存中的所有粒子索引
        if (cached.idxs[i] == GHOST_INDEX ||
            cached.idxs[i] == ATTRACTOR_INDEX) { // 如果是虚拟粒子或吸引子,则跳过
            continue; // ghost or attractor
        }
        cached.colors[i] = colorizer.evalColor(cached.idxs[i]); // 计算普通粒子的颜色
    }
}

static AutoPtr<IRenderContext> getContext(const RenderParams& params, Bitmap<Rgba>& bitmap) {//获取渲染上下文：渲染的前期准备接口类
    if (params.particles.doAntialiasing) {// 如果启用抗锯齿（FSAA）：消除画面中图形边缘的锯齿，使画面看起来更为平滑
        if (params.particles.smoothed) {//如果启用平滑化
            CubicSpline<2> kernel;// 使用CubicSpline二次核函数
            return makeAuto<SmoothedRenderContext>(bitmap, kernel);// 创建平滑过滤的渲染上下文 
        } else {
            return makeAuto<AntiAliasedRenderContext>(bitmap);// 否则创建基础抗锯齿渲染上下文
        }
    } else {//不需要FSAA
        if (params.background.a() == 1.f) {// 如果背景不透明
            return makeAuto<PreviewRenderContext<OverridePixelOp>>(bitmap);//创建预渲染上下文，进行OverridePixelOp操作：覆盖式绘制，用新颜色直接覆盖原像素颜色 
        } else {//背景透明
            return makeAuto<PreviewRenderContext<OverPixelOp>>(bitmap);//创建预渲染上下文，进行OverPixelOp操作：叠加式绘制，在原颜色上叠加新颜色
        }
    }
}

void ParticleRenderer::render(const RenderParams& params, Statistics& stats, IRenderOutput& output) const {//粒子渲染器的渲染逻辑
    MEASURE_SCOPE("ParticleRenderer::render");// 测量渲染总耗时

    Bitmap<Rgba> bitmap(params.camera->getSize());// 根据相机大小创建渲染目标位图
    AutoPtr<IRenderContext> context = getContext(params, bitmap);// 根据参数获取渲染上下文

    // fill with the background color
    context->fill(params.background);// 用背景颜色填充位图

    if (grid > 0.f) {// 如果指定了网格参数
        drawGrid(*context, *params.camera, grid);//绘制辅助网格
    }

    struct {
        Vector r;//位置向量
        Vector v;//着色器向量
        bool used = false;//使用标志，记录是否被选中
    } dir;//定义dir结构体用于临时存储单个粒子的位置和速度,从而实现对选中粒子的特殊处理

    context->setColor(Rgba::black(), ColorFlag::LINE);// 设置粒子默认轮廓颜色为黑色

    shouldContinue = true;//设置shouldContinue标志为true，表示继续渲染
    // draw particles 绘制粒子

    const bool reverseOrder = dot(cached.cameraDir, params.camera->getFrame().row(2)) < 0._f;//通过计算判断遍历粒子的顺序是否需要反转
    //前者cameraDir是相机视角的方向向量，后者是获取相机坐标系下的向上方向向量：getFrame() 获得相机的坐标系矩阵，row(2) 取矩阵的第3行,即获取了相机坐标系下的向上方向向量
    //计算二者的点积：如果点积小于0,表示两个向量方向相反,也就是相机视角翻转了180度，这时将 reverseOrder 设为 true。
    //在遍历粒子渲染时,如果 reverseOrder 为 true,则反转粒子遍历顺序，从最后一个粒子开始遍历，这样可以保证近端粒子绘制在远端粒子的上面，从而实现正确的粒子遮挡关系。

    for (Size k = 0; k < cached.positions.size(); ++k) {// 遍历所有粒子
        const Size i = reverseOrder ? cached.positions.size() - k - 1 : k; // 根据粒子遍历顺序（本质上是根据相机视角）,计算粒子索引
        if (!params.particles.renderGhosts && cached.idxs[i] == GHOST_INDEX) {// 如果是虚拟粒子且不渲染虚拟粒子,则跳过
            continue;
        }
        if (params.particles.selected && cached.idxs[i] == params.particles.selected.value()) {// 如果是选中的粒子
            // highlight the selected particle 高亮选中粒子
            context->setColor(Rgba::red(), ColorFlag::FILL);//设置填充颜色为红色
            context->setColor(Rgba::white(), ColorFlag::LINE);//设置轮廓颜色为白色

            if (!cached.vectors.empty()) {//如果缓存中的着色器向量数据不为空
                dir.used = true;//设置使用标志为true，代表被选中
                dir.v = cached.vectors[i];// 存储当前粒子的着色器向量
                dir.r = cached.positions[i];// 存储当前粒子的位置向量
            }
        } else {//如果不是选中粒子
            Rgba color = cached.colors[i];// 获取粒子的颜色
            if (params.particles.grayScale) {// 如果需要灰度化
                color = Rgba(color.intensity());// 转换为灰度颜色
            }
            context->setColor(color, ColorFlag::FILL | ColorFlag::LINE);// 设置粒子的填充颜色和轮廓颜色
            if (cached.idxs[i] == GHOST_INDEX || cached.idxs[i] == ATTRACTOR_INDEX) {//如果是虚拟粒子或吸引子
                context->setColor(Rgba::gray(0.7f), ColorFlag::LINE);// 设置浅灰色轮廓
            }
        }

        const Optional<ProjectedPoint> p = params.camera->project(cached.positions[i]);// 投影粒子到屏幕空间
        if (p) {//如果投影成功
            float radius;//粒子半径
            if (cached.idxs[i] != ATTRACTOR_INDEX) {// 如果不是吸引子
                radius = p->radius * params.particles.scale;// 计算粒子在屏幕空间的半径：按比例缩放
            } else {//是吸引子
                radius = p->radius;//不缩放
            }
            const float size = min<float>(radius, context->size().x);// 限制在画布最大尺寸以内
            context->drawCircle(p->coords, size);// 绘制粒子
        }
    }
    // after all particles are drawn, draw the velocity vector over 绘制完所有粒子后，绘制速度矢量
    if (dir.used) {//如果是被选中的粒子
        drawVector(*context, *params.camera, dir.r, dir.v, params.vectors.length);// 绘制选中粒子的速度向量
    }

    renderOverlay(*context, params, stats);// 渲染覆盖层

    // lastly black frame to draw on top of other stuff 最后绘制黑色边框,放在图片其他内容的上层
    const Pixel upper = bitmap.size() - Pixel(1, 1);// 计算图像大小,减去1个像素作为边框宽度
    context->setColor(Rgba::black(), ColorFlag::LINE);// 设置边框颜色为黑色
    //绘制矩形边框的四条边
    context->drawLine(Coords(0, 0), Coords(upper.x, 0));
    context->drawLine(Coords(upper.x, 0), Coords(upper));
    context->drawLine(Coords(upper), Coords(0, upper.y));
    context->drawLine(Coords(0, upper.y), Coords(0, 0));

    output.update(bitmap, context->getLabels(), true);// 更新输出
}

void ParticleRenderer::cancelRender() {//取消粒子渲染
    shouldContinue = false;//将shouldContinue设置为false来表示取消渲染
}

NAMESPACE_SPH_END
