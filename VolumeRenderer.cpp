#include "gui/renderers/VolumeRenderer.h"
#include "gui/Factory.h"
#include "gui/objects/Camera.h"
#include "gui/objects/Colorizer.h"
#include "gui/renderers/FrameBuffer.h"
#include "objects/finders/KdTree.h"
#include "objects/utility/OutputIterators.h"

NAMESPACE_SPH_BEGIN
//体积渲染：体积渲染是一组用于显示3D离散采样数据集（通常是3D标量场）的2D投影的技术。简而言之，就是将基于体素（体积的基本元素）的数据渲染为可见的2D图像。
//许多视觉效果，诸如烟雾、云、火焰，是很难用几何体表现的。这时候我们可以使用体积渲染技术，它通过计算光线在体积中的传播方式和相互作用，生成逼真的体积效果。
VolumeRenderer::VolumeRenderer(SharedPtr<IScheduler> scheduler, const GuiSettings& settings)//定义体积渲染器类
    : IRaytracer(scheduler, settings) {}//继承自光追类

VolumeRenderer::~VolumeRenderer() = default;//析构函数

const float MAX_DISTENTION = 50;//最大膨胀常量
const float MIN_NEIGHS = 8;//最小邻居常量

void VolumeRenderer::initialize(const Storage& storage, const IColorizer& colorizer, const ICamera& camera) {//初始化，接收Storage、IColorizer和ICamera作为参数
    cached.r = storage.getValue<Vector>(QuantityId::POSITION).clone();//从存储中获取位置向量,并克隆到cached.r。
    //这里缓存存储的是点的位置信息。
    this->setColorizer(colorizer);//设置着色器，用于将位置数据映射为RGB颜色。

    cached.distention.resize(cached.r.size());//存储粒子的膨胀系数。

    KdTree<KdNode> tree;
    //当我们场景中有大量物体的时候，每根射线都要和每个物体求交的消耗是非常昂贵的，我们需要对求交这个过程进行加速。
    //加速方法通常分为空间划分，物体划分。空间划分有八叉树、KDTree,BSPTree,物体划分有BVH、SAH等。
    tree.build(*scheduler, cached.r, FinderFlag::SKIP_RANK);//构建一个kd树（设置了跳过排序标志,可以加速构建过程）。
    //kd树：k-dimensional树的简称,构建KdTree时,会通过递归分割的方式,将点集划分为多个节点,每个节点下面再划分子节点,直到达到最大深度。
    //构建结果是一个包含层次节点的树结构, 通过这种空间划分,后续可以加速找到邻近点的查询。
    //用于后续的邻居点查找，可以加速查找位置周围的点, 从而优化渲染效率。
    Array<BvhSphere> spheres(cached.r.size());//用于存储每个点的包围球范围。
    spheres.reserve(cached.r.size());
    ThreadLocal<Array<NeighborRecord>> neighs(*scheduler);//使用 neighs 线程进行本地存储,为每个线程维护独立存储。
    parallelFor(*scheduler, neighs, 0, cached.r.size(), [&](const Size i, Array<NeighborRecord>& local) {//处理每个点
        const float initialRadius = cached.r[i][H];//计算初始半径
        float radius = initialRadius;
        while (radius < MAX_DISTENTION * initialRadius) {//从 initialRadius 开始扩大半径搜索邻居点。
            tree.findAll(cached.r[i], radius, local); //使用 KdTree 的 findAll 方法查找当前半径下的邻居点。
            if (local.size() >= MIN_NEIGHS) {//如果邻居点数 ≥ MIN_NEIGHS,则确定当前半径为最终半径。
                break;
            } else {
                radius *= 1.5f;//否则,将半径扩大1.5倍,继续搜索。
            }
        }

        BvhSphere s(cached.r[i], radius);//构建包围球数据，半径为确定的最终半径。
        s.userData = i;
        spheres[i] = s;// 将包围球写入 spheres 数组。

        cached.distention[i] = min(radius / initialRadius, MAX_DISTENTION);//计算距离膨胀系数：取最终半径与初始半径的比值和最大膨胀常量中的最小值。
    });

    //基于质量和材质信息估算出每个点的合理参考半径
    ArrayView<const Float> m = storage.getValue<Float>(QuantityId::MASS);// 获取点质量数据m。
    cached.referenceRadii.resize(cached.r.size());//调整参考半径数组的大小与点的数量相同。
    if (storage.getMaterialCnt() > 0) {
        for (Size matId = 0; matId < storage.getMaterialCnt(); ++matId) {//对于每个材质
            MaterialView mat = storage.getMaterial(matId);
            const Float rho = mat->getParams().has(BodySettingsId::DENSITY)//获取材质密度rho。
                                  ? mat->getParam<Float>(BodySettingsId::DENSITY)
                                  : 1000._f;
            for (Size i : mat.sequence()) {//对属于该材质的点i:
                const Float volume = m[i] / rho;//计算点的体积 = 质量 / 密度
                cached.referenceRadii[i] = root<3>(3._f * volume / (4._f * PI));//通过体积计算参考半径
            }
        }
    } else {//如果没有材质信息
        // guess the dentity
        const Float rho = 1000._f;//使用默认密度1000，估算参考半径
        for (Size i = 0; i < m.size(); ++i) {
            const Float volume = m[i] / rho;
            cached.referenceRadii[i] = root<3>(3._f * volume / (4._f * PI));
        }
    }

    cached.textures.clear();//清空纹理数组
    for (Size i = 0; i < storage.getAttractorCnt(); ++i) {//遍历所有吸引子attractor
        const Attractor& a = storage.getAttractors()[i];//从storage获取单个attractor
        const bool visible = a.settings.getOr(AttractorSettingsId::VISIBLE, true);//读取attractor的可见性参数
        const Float albedo = a.settings.getOr(AttractorSettingsId::ALBEDO, 1._f);//读取attractor的反射率参数
        cached.attractors.push(AttractorData{ a.mass, a.position, a.radius, visible, albedo });//将吸引子参数存入cached.attractors数组。

        String texturePath = a.settings.getOr<String>(AttractorSettingsId::VISUALIZATION_TEXTURE, "");//获取attractor的纹理路径
        if (!texturePath.empty()) {//如果有纹理
            if (cached.textureCache.contains(texturePath)) {//如果缓存中已存在,直接push进cached.textures数组
                cached.textures.push(cached.textureCache[texturePath]);
            } else {//否则加载纹理,放入缓存并push。
                SharedPtr<Texture> texture =
                    makeShared<Texture>(Path(texturePath), TextureFiltering::BILINEAR);
                cached.textureCache.insert(texturePath, texture);
                cached.textures.push(texture);
            }
        } else {// 如果无纹理,push一个nullptr
            cached.textures.push(nullptr);
        }

        BvhSphere sphere(a.position, a.radius);//创建attractor的包围球,放入spheres数组。
        sphere.userData = cached.r.size() + i;
        spheres.push(sphere);
    }

    bvh.build(std::move(spheres));//.使用包围球数据构建加速结构bvh。
    //层次包围体BVH（Bounding Volume Hierachies）：BV（Bounding Volume）是包含一组物体的空间体，它比其所包含的几何物体形状要简单的多，所以对它进行碰撞检测速度比直接与物体本身求相交更快。
    //BVH是一种以物体BV为基础进行划分的结构。本质上是一棵树。它由根节点、内部节点和叶子节点组成。其中叶子节点存放物体，每个非叶子节点都有包围体，父节点可以把子节点包围起来。
    
    //计算场景最大距离
    cached.maxDistance = 0;
    for (const Attractor& a : storage.getAttractors()) {//遍历所有attractor
        const Float dist = getLength(a.position - camera.getPosition());//对每个attractor,计算它到相机的距离dist。
        cached.maxDistance = max(cached.maxDistance, 2 * dist);//将dist的2倍与当前的cached.maxDistance对比,保留较大值作为新的cached.maxDistance。
    }

    for (ThreadData& data : threadData) {//初始化线程数据threadData
        data.data = RayData{};
    }

    shouldContinue = true;//设置渲染状态shouldContinue为true，确保渲染可以开始
}

bool VolumeRenderer::isInitialized() const {//检查缓存数据是否为空，以此判断是否进行了初始化。
    return !cached.r.empty();
}

void VolumeRenderer::setColorizer(const IColorizer& colorizer) {//设置着色器
    cached.colors.resize(cached.r.size());//调整颜色数组大小
    for (Size i = 0; i < cached.r.size(); ++i) {//调用evalColor计算颜色
        cached.colors[i] = colorizer.evalColor(i);
    }
}

Rgba VolumeRenderer::shade(const RenderParams& params, const CameraRay& cameraRay, ThreadData& data) const {//阴影着色
    const Vector primaryDir = getNormalized(cameraRay.target - cameraRay.origin);//计算摄像机射线的主方向primaryDir:从origin指向target
    const Ray primaryRay(cameraRay.origin, primaryDir);//使用主方向构建主射线primaryRay

    RayData& rayData(data.data);//线程数据：包含了渲染数据
    LensingEffect::Segments& segments = rayData.segments;//射线段数据
    Array<CurvedRayIntersectionInfo>& intersections = rayData.intersections;//交点数据

    /// 调用LensingEffect类,基于透镜效应获取所有的交点。
    LensingEffect lensing(cached.attractors,
        params.relativity.lensingMagnitude,
        0.1_f * cached.maxDistance,
        cached.maxDistance,
        params.volume.absorption > 0.f);
    //输入:attractors、lensing参数、场景最大距离、吸收系数，构造lensing对象
    //输出: 射线段segments和交点intersections
    const Ray lastRay = lensing.getAllIntersections(bvh, primaryRay, segments, intersections);
    //调用LensingEffect的getAllIntersections方法,基于主射线primaryRay计算重力透镜效应的所有交点和射线段。存储到segments和intersections中。
    //得到最终射线lastRay
    Rgba result = this->getEnviroColor(CameraRay{ lastRay.origin(), lastRay.origin() + lastRay.direction() });//基于lastRay计算环境光颜色

    for (const CurvedRayIntersectionInfo& is : reverse(intersections)) {// 反向遍历所有光线交点
        const BvhSphere* s = static_cast<const BvhSphere*>(is.object);// 将交点对象转换为球体指针
        const Size i = s->userData;// 获取球体的userData作为粒子索引
        const Vector hit = is.segment->origin() + is.segment->direction() * is.t; // 计算碰撞点坐标 = 光线原点 + 方向 * 交点距离
        const Vector center = s->getCenter();// 获取球心坐标
        const Vector toCenter = getNormalized(center - hit);// 计算碰撞点指向球心方向的向量，经过归一化作为单位向量
        
        //处理光线碰到attractor的情况
        //吸引子是一种特殊的粒子，它会产生吸引域，对一般粒子产生影响作用。
        if (i >= cached.r.size()) {//// 如果索引超过粒子数目,则为attractor
            // attractor, a solid object -> erase the emission accumulated so far if visible  吸引子是一类固体，需要通过提前着色，清除目前累积的可见的发光量
            const Size idx = i - cached.r.size();// 吸引子索引 = 交点索引 - 粒子数目
            if (cached.attractors[idx].visible) {//如果可见
                result = this->getAttractorColor(params, idx, hit);//调用getAttractorColor计算attractor在碰撞点hit的颜色
            }
            continue;
        }
        // 否则是粒子,进行标准着色计算
        const float cosPhi = abs(dot(toCenter, is.segment->direction()));// 计算碰撞点到球心向量与光线方向的余弦值
        const float distention = cached.distention[i]; //获取第i个粒子的膨胀系数
        // smoothing length should not have effect on the total emission 平滑长度不应对总发光量产生影响
        const float radiiFactor = cached.referenceRadii[i] / cached.r[i][H];//计算粒子半径因子radiiFactor,用参考半径除以当前半径,以消除光滑长度的影响。
        const float secant = 2._f * getLength(center - hit) * cosPhi * radiiFactor;//计算交点到球心的割线长度secant,用于后续计算发光强度
        // make dilated particles absorb more 使膨胀的粒子吸收更多光线
        result = result * exp(-params.volume.absorption * secant * pow<3>(cosPhi/distention));// 使用 secant 和 cosPhi 计算吸收量,吸收随膨胀程度的增大而增强
        // 3th power of cosPhi to give more weight to the sphere center,
        // divide by distention^3; distention should not affect the total emission  乘以cospi的3次方是为了使中心的影响权重更大；除以膨胀系数的三次方是为了使膨胀不影响总发光量。
        const float magnitude = params.volume.emission * pow<3>(cosPhi / distention) * secant;// 计算发光量，中心权重高,不受膨胀影响
        result += cached.colors[i] * magnitude;//将发光量乘以粒子颜色得到该粒子的发光效果，并将其添加到结果中。
        result.a() += magnitude;//将发光量添加到alpha通道：阿尔法通道(Alpha Channel)是指一张图片的透明和半透明度
    }
    result.a() = min(result.a(), 1.f);// 将最终alpha限制在1以内
    return result;// 返回结果
}

Rgba VolumeRenderer::getAttractorColor(const RenderParams& params,//吸引子着色
    const Size index,
    const Vector& hit) const {
    const AttractorData& a = cached.attractors[index];//// 获取当前attractor数据
    Rgba diffuse = Rgba::gray(a.albedo);// 初始化漫反射颜色为灰色
    const SharedPtr<Texture>& texture = cached.textures[index];//// 获取attractor纹理
    if (texture) {//如果有纹理
        const Vector r0 = hit - a.position;//计算碰撞点在attractor坐标系下的坐标向量
        SphericalCoords spherical = cartensianToSpherical(r0);//将坐标向量转换为球坐标
        Vector uvw = Vector(0.5_f - spherical.phi / (2._f * PI), spherical.theta / PI, 0._f);// 将球坐标转换为uvw纹理坐标
        diffuse = texture->eval(uvw) * a.albedo;//采样纹理，与albedo漫反射系数相乘作为最终漫反射颜色
    }

    const Vector n = getNormalized(a.position - hit);//计算碰撞点指向attractor方向，经过归一化作为单位向量
    const Float cosPhi = dot(n, params.lighting.dirToSun);//计算n向量与光照方向夹角的余弦值cosPhi
    if (cosPhi <= 0._f) {// 如果余弦<=0,光线与n向量反向或垂直 
        // not illuminated -> just ambient light  不会被光线照射到,仅返回环境光
        return diffuse * params.lighting.ambientLight;//返回环境光颜色
    }
        return diffuse * float(PI * cosPhi * params.lighting.sunLight + params.lighting.ambientLight);//否则返回: 漫反射颜色*(环境光照 + 太阳光照 * 余弦)
   }

NAMESPACE_SPH_END
