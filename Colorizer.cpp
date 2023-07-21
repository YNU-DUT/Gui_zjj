#include "gui/objects/Colorizer.h"
#include "physics/Eos.h"
#include "sph/Materials.h"
#include "timestepping/TimeStepCriterion.h"

NAMESPACE_SPH_BEGIN

DirectionColorizer::DirectionColorizer(const Vector& axis, const Palette& palette)//方向着色器构造函数
    : palette(palette)
    , axis(axis) {
    SPH_ASSERT(almostEqual(getLength(axis), 1._f));// 判断axis是否为单位向量
    // compute 2 perpendicular directions 计算与axis垂直的两个单位向量dir1和dir2
    Vector ref;
    if (almostEqual(axis, Vector(0._f, 0._f, 1._f)) || almostEqual(axis, Vector(0._f, 0._f, -1._f))) {//如果axis的方向为z轴
        ref = Vector(0._f, 1._f, 0._f);//定义ref为y轴方向上的单位向量
    } else {
        ref = Vector(0._f, 0._f, 1._f);//否则定义ref为z轴方向上的单位向量
    }
    dir1 = getNormalized(cross(axis, ref));// 通过axis和ref叉乘计算dir1
    dir2 = cross(axis, dir1); // 通过axis和dir1叉乘计算dir2
    SPH_ASSERT(almostEqual(getLength(dir2), 1._f));// 判断dir2是否为单位向量
}

Optional<float> DirectionColorizer::evalScalar(const Size idx) const {//根据索引计算向量的标量值
    SPH_ASSERT(this->isInitialized());// 检查是否初始化
    const Vector projected = values[idx] - dot(values[idx], axis) * axis;// 将values投影到axis的法向量平面上
    //计算values在axis方向上的投影长度，生成一个在axis方向上的投影向量，从values向量中抵消这个在axis方向上的投影向量，从而消除了values在axis方向的分量,使其投影到与axis垂直的平面上
    const Float x = dot(projected, dir1);//投影向量projected在dir1方向上的投影长度
    const Float y = dot(projected - x * dir1, dir2);//在dir2方向上的投影长度
    //计算出了values向量在与axis垂直平面上即参考系(dir1, dir2)上的坐标(x, y)
    return float(PI + atan2(y, x));//将(x,y)坐标转换为极角，将二维坐标映射到一维空间
    // 极坐标下, 角度和颜色之间的映射是均匀连续的,便于后续进行颜色映射
}

static thread_local Array<NeighborRecord> neighs;

SummedDensityColorizer::SummedDensityColorizer(const RunSettings& settings, Palette palette)//密度贡献着色器构造函数
    : palette(std::move(palette)) {
    finder = Factory::getFinder(settings);//创建密度点查找器finder
    kernel = Factory::getKernel<3>(settings);//创建3维核函数求解器kernel
}

void SummedDensityColorizer::initialize(const Storage& storage, const RefEnum ref) {//密度贡献着色器初始化
    //从storage中取出质量和位置数据到m和r
    m = makeArrayRef(storage.getValue<Float>(QuantityId::MASS), ref);
    r = makeArrayRef(storage.getValue<Vector>(QuantityId::POSITION), ref);

    finder->build(SEQUENTIAL, r);//用位置数据r构建密度点查找器finder
}

float SummedDensityColorizer::sum(const Size idx) const {//累加密度贡献值，计算密度和
    finder->findAll(idx, r[idx][H] * kernel.radius(), neighs);// 根据索引值idx在半径范围内找到附近的所有相邻点
    Float rho = 0._f;// 初始化密度和为0
    for (const auto& n : neighs) {// 遍历所有相邻点
        rho += m[n.index] * kernel.value(r[idx] - r[n.index], r[idx][H]);// 累加每个相邻点的密度贡献值
    }
    return float(rho);
}

bool TemperatureColorizer::hasData(const Storage& storage) const {// 检查存储数据是否满足温度着色的需求
    if (!storage.has(QuantityId::DENSITY) || !storage.has(QuantityId::ENERGY)) {// // 检查是否有密度和能量数据
        return false;
    }
    if (storage.getMaterialCnt() == 0) {// 检查是否有材料数据
        return false;
    }
    for (Size matId = 0; matId < storage.getMaterialCnt(); ++matId) { // 检查材料是否都是EosMaterial类型
        const IMaterial& material = storage.getMaterial(matId);
        if (!dynamic_cast<const EosMaterial*>(&material)) {
            return false;
        }
    }
    return true;
}

void TemperatureColorizer::initialize(const Storage& storage, const RefEnum UNUSED(ref)) {//温度着色器初始化
    // 获取密度和能量数据
    ArrayView<const Float> rho = storage.getValue<Float>(QuantityId::DENSITY);
    ArrayView<const Float> u = storage.getValue<Float>(QuantityId::ENERGY);
    Array<Float> T(rho.size());// 获取密度和能量数据
    for (Size matId = 0; matId < storage.getMaterialCnt(); ++matId) {// 对每个材料
        MaterialView view = storage.getMaterial(matId); // 获取材料视图
        EosMaterial& material = dynamic_cast<EosMaterial&>(view.material()); // 获取材料对象
        for (Size i : view.sequence()) {// 对材料中的每个点
            T[i] = material.getEos().getTemperature(rho[i], u[i]);// 计算温度
        }
    }
    values = makeArrayRef(std::move(T));//将要着色的值存入温度数组
}

void CorotatingVelocityColorizer::initialize(const Storage& storage, const RefEnum ref) {//旋转速度着色器初始化
    // 获取位置、速度和材料id
    r = makeArrayRef(storage.getValue<Vector>(QuantityId::POSITION), ref);
    v = makeArrayRef(storage.getDt<Vector>(QuantityId::POSITION), ref);
    matIds = makeArrayRef(storage.getValue<Size>(QuantityId::MATERIAL_ID), ref);

    ArrayView<const Float> m = storage.getValue<Float>(QuantityId::MASS);// 获取质量
    // 将位置数据和速度数据转换为数组视图
    ArrayView<const Vector> rv = r;
    ArrayView<const Vector> vv = v;
    data.resize(storage.getMaterialCnt());// 初始化材料数据数组

    for (Size i = 0; i < data.size(); ++i) {// 为每个材料计算质心和角速度
        MaterialView mat = storage.getMaterial(i);// 获取材料视图
        // 计算该材料粒子的索引范围
        const Size from = *mat.sequence().begin();
        const Size to = *mat.sequence().end();
        const Size size = to - from;// 计算该材料的粒子数量
        //  使用getCenterOfMass函数计算该材料的质心
        data[i].center = Post::getCenterOfMass(m.subset(from, size), rv.subset(from, size));
        //  使用getAngularFrequency函数计算该材料的角速度
        data[i].omega =
            Post::getAngularFrequency(m.subset(from, size), rv.subset(from, size), vv.subset(from, size));
    }
}

bool DamageActivationColorizer::hasData(const Storage& storage) const {// 检查存储数据是否满足损伤激活着色的需求
    // 检查是否存在以下数据:
    // 1. 应力张量的偏微分(Deviatoric stress)
    // 2. 压强(Pressure)
    // 3. 最小主应变(Eps min)
    // 4. 损伤值(Damage)
    return storage.has(QuantityId::DEVIATORIC_STRESS) && storage.has(QuantityId::PRESSURE) &&
           storage.has(QuantityId::EPS_MIN) && storage.has(QuantityId::DAMAGE);
}

void DamageActivationColorizer::initialize(const Storage& storage, const RefEnum UNUSED(ref)) {//损伤激活着色器初始化
    // 获取应力张量偏微分、压强、最小主应变、损伤值
    ArrayView<const TracelessTensor> s = storage.getValue<TracelessTensor>(QuantityId::DEVIATORIC_STRESS);
    ArrayView<const Float> p = storage.getValue<Float>(QuantityId::PRESSURE);
    ArrayView<const Float> eps_min = storage.getValue<Float>(QuantityId::EPS_MIN);
    ArrayView<const Float> damage = storage.getValue<Float>(QuantityId::DAMAGE);

    ratio.resize(p.size());// 初始化比率数组
    /// \todo taken from ScalarGradyKippDamage, could be deduplicated 取自ScalarGradyKippDamage，可以去重
    for (Size matId = 0; matId < storage.getMaterialCnt(); ++matId) { // 遍历每个材料
        MaterialView mat = storage.getMaterial(matId);// 获取材料视图
        const Float young = mat->getParam<Float>(BodySettingsId::YOUNG_MODULUS); // 获取杨氏模量

        /// \todo parallelize 并行化
        for (Size i : mat.sequence()) {// 遍历材料中的每个粒子
            const SymmetricTensor sigma = SymmetricTensor(s[i]) - p[i] * SymmetricTensor::identity();// 计算粒子的应力张量
            // 求解应力张量的特征值
            Float sig1, sig2, sig3;
            tie(sig1, sig2, sig3) = findEigenvalues(sigma);
            const Float sigMax = max(sig1, sig2, sig3);// 最大特征值
            const Float young_red = max((1._f - pow<3>(damage[i])) * young, 1.e-20_f);// 计算减小的杨氏模量
            const Float strain = sigMax / young_red;// 计算最大应变
            ratio[i] = float(strain / eps_min[i]);// 计算比率
        }
    }
}

BeautyColorizer::BeautyColorizer() {//美化着色器构造函数
    palette = Palette({ { 0.f, Rgba(0.5f, 0.5f, 0.5) }, { 1.f, Rgba(1.f, 1.f, 0.6f) } },// 初始化调色板
        Interval(u_0, u_yellow),
        PaletteScale::LOGARITHMIC);//对数尺度
    /// \todo maybe not fixed? 能不能不固定
    // 在调色板上添加两个固定颜色点
    palette.addFixedPoint(u_glow, Rgba(0.5f, 0.5f, 0.5f));//发光
    palette.addFixedPoint(u_red, Rgba(0.8f, 0.f, 0.f));//红色
    f_glow = (log10(u_glow) - log10(u_0)) / (log10(u_yellow) - log10(u_0));//计算发光区间的映射因子：发光区间为[u_0, u_glow]，通过对数尺度下的区间长度比例计算映射因子
    //计算映射因子是为了对发光区间做线性映射, 使其在最终颜色结果中有明显体现。
}

BoundaryColorizer::BoundaryColorizer(const Detection detection, const Float threshold)//边界着色器构造函数
    : detection(detection) {
    if (detection == Detection::NEIGBOUR_THRESHOLD) {// 相邻数阈值检测:具有较少邻居的粒子被认为是边界粒子。增加阈值会导致边界粒子增多。
        neighbors.threshold = Size(threshold);//设定相邻数阈值
    } else {// 法线阈值检测:边界由近似曲面法线的相对位置向量确定。开销较高，但不敏感地依赖于邻居粒子的数量。增加阈值会导致边界粒子减少。
        normals.threshold = threshold;//设定法线阈值
    }
}

bool BoundaryColorizer::hasData(const Storage& storage) const {// 检查存储数据是否满足边界着色的需求
    if (detection == Detection::NORMAL_BASED) {// 如果检测方式为法线
        return storage.has(QuantityId::SURFACE_NORMAL);// 检查是否存在表面法线数据
    } else {// 如果检测方式为相邻数
        return storage.has(QuantityId::NEIGHBOR_CNT);// 检查是否存在相邻数数据
    }
}

void BoundaryColorizer::initialize(const Storage& storage, const RefEnum ref) {//边界着色器初始化
    // 根据检测方式获取对应数据
    if (detection == Detection::NORMAL_BASED) {// 如果检测方式为法线
        normals.values = makeArrayRef(storage.getValue<Vector>(QuantityId::SURFACE_NORMAL), ref);// 获取表面法线数据
    } else {// 如果检测方式为相邻数
        neighbors.values = makeArrayRef(storage.getValue<Size>(QuantityId::NEIGHBOR_CNT), ref);// 获取相邻数数据
    }
}

bool BoundaryColorizer::isInitialized() const {//边界着色器初始化判断
    return (detection == Detection::NORMAL_BASED && !normals.values.empty()) ||//// 根据检测方式检查对应数据
           (detection == Detection::NEIGBOUR_THRESHOLD && !neighbors.values.empty());
}

Rgba BoundaryColorizer::evalColor(const Size idx) const {//边界着色器的分类着色
    if (isBoundary(idx)) {// 调用isBoundary函数判断索引idx对应的粒子是否在边界
    //isBoundary函数根据构造函数设置的边界检测方式,使用阈值判断idx对应的粒子是否在边界
        return Rgba::red();//如果是边界粒子,着色为红色
    } else {
        return Rgba::gray();//如果不是边界粒子,着色为灰色
    }
}

bool BoundaryColorizer::isBoundary(const Size idx) const {//边界判断
    switch (detection) {
    case Detection::NEIGBOUR_THRESHOLD:
        SPH_ASSERT(!neighbors.values.empty());//确保有数据
        return neighbors.values[idx] < neighbors.threshold;//检查相邻数是否低于阈值
    case Detection::NORMAL_BASED:
        SPH_ASSERT(!normals.values.empty());
        return getLength(normals.values[idx]) > normals.threshold;//检查法线长度是否大于阈值
    default:
        NOT_IMPLEMENTED;
    }
}


/// \todo possibly move elsewhere 可能会搬到其他地方
static uint64_t getHash(const uint64_t value, const Size seed) {//生成一个64位无符号整数的哈希值
    // https://stackoverflow.com/questions/8317508/hash-function-for-a-string
    constexpr int A = 54059;
    constexpr int B = 76963;
    constexpr int FIRST = 37;

    uint64_t hash = FIRST + seed;
    StaticArray<uint8_t, sizeof(uint64_t)> data;
    std::memcpy(&data[0], &value, data.size());
    for (Size i = 0; i < sizeof(uint64_t); ++i) {
        hash = (hash * A) ^ (data[i] * B);
    }
    return hash;
}

static Rgba getRandomizedColor(const Size idx, const Size seed = 0) {//利用hash高质量的随机性，根据输入的索引和种子生成一个随机颜色值
    const uint64_t hash = getHash(idx, seed); //使用getHash函数根据idx和seed生成一个64位无符号整数的哈希值hash
    //从hash中提取出红色、绿色和蓝色的随机值
    //r是取hash低16位
    //g是取hash中16 - 32位
    //b是取hash高32位
    const uint8_t r = (hash & 0x00000000FFFF);
    const uint8_t g = (hash & 0x0000FFFF0000) >> 16;
    const uint8_t b = (hash & 0xFFFF00000000) >> 32;
    //将取出的随机值映射到0 - 255的范围内,转成相应的红绿蓝分量。使用红绿蓝分量构造出一个Rgba颜色值并返回。
    return Rgba(r / 255.f, g / 255.f, b / 255.f);
}

template <typename TDerived>//模板函数
Rgba IdColorizerTemplate<TDerived>::evalColor(const Size idx) const {//IdColorizer的着色
    const Optional<Size> id = static_cast<const TDerived*>(this)->evalId(idx);//调用evalId()函数根据索引值idx获取对应的id值
    if (!id) {
        return Rgba::gray();//如果id为空,直接返回灰色
    }
    const Rgba color = getRandomizedColor(id.value(), seed);// 若id存在,则调用getRandomizedColor根据id和种子seed生成随机颜色
    return color;
}

template <typename TDerived>
Optional<Particle> IdColorizerTemplate<TDerived>::getParticle(const Size idx) const {//获取粒子
    Particle particle(idx);//创建一个粒子对象
    const Optional<Size> id = static_cast<const TDerived*>(this)->evalId(idx);//获取id
    if (id) {
        particle.addValue(QuantityId::FLAG, id.value());//若id存在，则设置粒子的FLAG属性设置为idx值
    }
    return particle;
}

//模板类
template class IdColorizerTemplate<ParticleIdColorizer>;//粒子ID着色器
template class IdColorizerTemplate<ComponentIdColorizer>;//组件ID着色器
template class IdColorizerTemplate<AggregateIdColorizer>;//集群ID着色器
template class IdColorizerTemplate<IndexColorizer>;//索引着色器

void ParticleIdColorizer::initialize(const Storage& storage, const RefEnum ref) {//粒子ID着色器初始化
    if (storage.has(QuantityId::PERSISTENT_INDEX)) {//若存储中存在持续索引
        persistentIdxs = makeArrayRef(storage.getValue<Size>(QuantityId::PERSISTENT_INDEX), ref);//从存储中获取持续索引,并存储到类中
    }
}

Optional<Particle> ParticleIdColorizer::getParticle(const Size idx) const {//粒子ID着色器获取粒子
    Particle particle(idx);
    particle.addValue(QuantityId::FLAG, idx);//设置粒子的FLAG属性设置为idx值
    if (!persistentIdxs.empty() && idx < persistentIdxs.size()) {//若持久索引数组不为空且索引有效
        particle.addValue(QuantityId::PERSISTENT_INDEX, persistentIdxs[idx]);//取出对应持久索引,添加到粒子的对应属性中
    }
    return particle;
}


ComponentIdColorizer::ComponentIdColorizer(const GuiSettings& gui,//// 组件ID着色器构造函数
    const Flags<Post::ComponentFlag> connectivity,
    const Optional<Size> highlightIdx)
    : IdColorizerTemplate<ComponentIdColorizer>(gui)
    , connectivity(connectivity)//存储组件连接性信息
    , highlightIdx(highlightIdx) {}//存储要高亮的组件索引

void ComponentIdColorizer::setHighlightIdx(const Optional<Size> newHighlightIdx) {// 设置高亮组件索引
    if (newHighlightIdx) {
        highlightIdx = min(newHighlightIdx.value(), components.size() - 1);//如果新索引有效,取它和组件数目的最小值作为高亮索引
    } else {
        highlightIdx = NOTHING;//如果新索引为空,则设置高亮索引为无效值
    }
}

Rgba ComponentIdColorizer::evalColor(const Size idx) const {//组件ID着色器着色
    if (highlightIdx) {// 如果存在高亮索引
        if (highlightIdx.value() == components[idx]) {//若当前组件为高亮组件
            return Rgba(1.f, 0.65f, 0.f);//高亮颜色
        } else {
            return Rgba::gray(0.3f);//否则设为背景颜色
        }
    } else { // 没有高亮索引 
        const Size id = this->evalId(idx).value();// 获取组件id
        return getRandomizedColor(compIdxs[id], seed);// 根据id选择随机颜色
    }
}

Optional<Particle> ComponentIdColorizer::getParticle(const Size idx) const {// 组件ID着色器获取粒子
    Particle particle(idx);// 创建粒子对象
    const Optional<Size> id = this->evalId(idx); // 获取组件id
    particle.addValue(QuantityId::FLAG, id.value());//设置粒子的FLAG属性
    // 找到属于同一组件的所有粒子索引
    Array<Size> indices;
    for (Size i = 0; i < r.size(); ++i) {
        if (components[i] == id.value()) {
            indices.push(i);
        }
    }
    if (indices.size() > 1) { // 如果组件含多个粒子,计算组件的角速度
        const Vector omega = Post::getAngularFrequency(m, r, v, indices);
        particle.addValue(QuantityId::ANGULAR_FREQUENCY, getLength(omega));//设置粒子的ANGULAR_FREQUENCY属性
    }
    return particle;
}

bool ComponentIdColorizer::hasData(const Storage& storage) const {// 检查存储数据是否满足组件ID着色的需求
    return hasVelocity(storage);// 检查存储中是否存在速度数据
}

void ComponentIdColorizer::initialize(const Storage& storage, const RefEnum ref) {//组件ID着色器初始化
    const Array<Vector>& current = storage.getValue<Vector>(QuantityId::POSITION); //获取当前位置数据
    if (current == cached.r) {//和缓存位置数据比较：如果相同,跳过初始化;如果不同,则进行初始化
        // optimization, very poorly done 优化，做的很差
        return;
    }

    m = makeArrayRef(storage.getValue<Float>(QuantityId::MASS), ref);//获取质量数据
    r = makeArrayRef(storage.getValue<Vector>(QuantityId::POSITION), ref);//获取位置数据
    v = makeArrayRef(storage.getDt<Vector>(QuantityId::POSITION), ref);//获取速度数据

    cached.r = current.clone();//将当前位置数据current复制存储到缓存cached.r中

    const Size numComponents = Post::findComponents(storage, 2.5_f, connectivity, components);// 使用findComponents函数查找粒子组件

    // sort by the smallest index in each component 按每个组件中最小的索引排序
    //初始化compIdxs数组为numComponents大小, 填充为-1
    compIdxs.resize(numComponents);
    compIdxs.fill(Size(-1));
    for (Size i = 0; i < components.size(); ++i) {//遍历components数组，记录组件的最小粒子索引
        const Size ci = components[i];
        compIdxs[ci] = min(compIdxs[ci], i);
    }
}

String ComponentIdColorizer::name() const {//根据connectivity的设置返回不同的名称,以反映组件的划分方式
    if (connectivity.has(Post::ComponentFlag::ESCAPE_VELOCITY)) {//检查连接参数中是否设置了ESCAPE_VELOCITY（按逃逸速度划分）属性
        return "Bound component ID";//返回"Bound component ID"
    } else if (connectivity.has(Post::ComponentFlag::SEPARATE_BY_FLAG)) {//检查连接参数中是否设置了SEPARATE_BY_FLAG（按FLAG划分）属性
        return "Component ID (flag)";//返回"Component ID (flag)"
    } else {//如果都未设置
        return "Component ID";//返回默认"Component ID"
    }
}

void MaterialColorizer::initialize(const Storage& storage, const RefEnum ref) {//材料着色器初始化
    IndexColorizer::initialize(storage, ref);// 调用父类的initialize做基础初始化

    const Size matCnt = storage.getMaterialCnt();// 获取材料总数
    //初始化存储EOS和Rheology类型的数组
    //EOS：Equation of state for this material 状态方程
    //Rheology：Model of stress reducing used within the rheological model 流变模型中使用的应力减小模型
    eosNames.resize(matCnt);
    rheoNames.resize(matCnt);
    for (Size matId = 0; matId < matCnt; ++matId) {//遍历每个材料
        const IMaterial& mat = storage.getMaterial(matId);//获取材料信息
        //获取EOS和Rheology类型，并存入对应数组
        const EosEnum eos = mat.getParam<EosEnum>(BodySettingsId::EOS);
        const YieldingEnum yield = mat.getParam<YieldingEnum>(BodySettingsId::RHEOLOGY_YIELDING);
        eosNames[matId] = EnumMap::toString(eos);
        rheoNames[matId] = EnumMap::toString(yield);
    }
}

Optional<Particle> MaterialColorizer::getParticle(const Size idx) const {//材料着色器获取粒子
    Particle particle(idx);
    const Optional<Size> id = IndexColorizer::evalId(idx);// // 获取材料id
    if (id) {//id存在
        particle.addValue(QuantityId::MATERIAL_ID, id.value());//设置粒子的MATERIAL_ID属性
        //设置粒子的EOS和RHEOLOGY_YIELDING属性
        particle.addParameter(BodySettingsId::EOS, eosNames[id.value()]);
        particle.addParameter(BodySettingsId::RHEOLOGY_YIELDING, rheoNames[id.value()]);
    }
    return particle;
}

void TimeStepColorizer::initialize(const Storage& storage, const RefEnum ref) {//时间步长着色器初始化
    TypedColorizer<Float>::initialize(storage, ref);//调用父类的initialize做基础初始化
    critIds = makeArrayRef(storage.getValue<Size>(QuantityId::TIME_STEP_CRITERION), ref);//获取时间步长判定标准
}

Optional<Particle> TimeStepColorizer::getParticle(const Size idx) const {//时间步长着色器获取粒子
    Particle particle(idx);
    particle.addValue(QuantityId::TIME_STEP, values[idx]);// 设置粒子的时间步长属性
    particle.addValue(QuantityId::TIME_STEP_CRITERION, toString(CriterionId(critIds[idx])));// 设置粒子的时间步长判定标准属性
    return particle;
}

NAMESPACE_SPH_END
