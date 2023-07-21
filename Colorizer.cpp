#include "gui/objects/Colorizer.h"
#include "physics/Eos.h"
#include "sph/Materials.h"
#include "timestepping/TimeStepCriterion.h"

NAMESPACE_SPH_BEGIN

DirectionColorizer::DirectionColorizer(const Vector& axis, const Palette& palette)//������ɫ�����캯��
    : palette(palette)
    , axis(axis) {
    SPH_ASSERT(almostEqual(getLength(axis), 1._f));// �ж�axis�Ƿ�Ϊ��λ����
    // compute 2 perpendicular directions ������axis��ֱ��������λ����dir1��dir2
    Vector ref;
    if (almostEqual(axis, Vector(0._f, 0._f, 1._f)) || almostEqual(axis, Vector(0._f, 0._f, -1._f))) {//���axis�ķ���Ϊz��
        ref = Vector(0._f, 1._f, 0._f);//����refΪy�᷽���ϵĵ�λ����
    } else {
        ref = Vector(0._f, 0._f, 1._f);//������refΪz�᷽���ϵĵ�λ����
    }
    dir1 = getNormalized(cross(axis, ref));// ͨ��axis��ref��˼���dir1
    dir2 = cross(axis, dir1); // ͨ��axis��dir1��˼���dir2
    SPH_ASSERT(almostEqual(getLength(dir2), 1._f));// �ж�dir2�Ƿ�Ϊ��λ����
}

Optional<float> DirectionColorizer::evalScalar(const Size idx) const {//�����������������ı���ֵ
    SPH_ASSERT(this->isInitialized());// ����Ƿ��ʼ��
    const Vector projected = values[idx] - dot(values[idx], axis) * axis;// ��valuesͶӰ��axis�ķ�����ƽ����
    //����values��axis�����ϵ�ͶӰ���ȣ�����һ����axis�����ϵ�ͶӰ��������values�����е��������axis�����ϵ�ͶӰ�������Ӷ�������values��axis����ķ���,ʹ��ͶӰ����axis��ֱ��ƽ����
    const Float x = dot(projected, dir1);//ͶӰ����projected��dir1�����ϵ�ͶӰ����
    const Float y = dot(projected - x * dir1, dir2);//��dir2�����ϵ�ͶӰ����
    //�������values��������axis��ֱƽ���ϼ��ο�ϵ(dir1, dir2)�ϵ�����(x, y)
    return float(PI + atan2(y, x));//��(x,y)����ת��Ϊ���ǣ�����ά����ӳ�䵽һά�ռ�
    // ��������, �ǶȺ���ɫ֮���ӳ���Ǿ���������,���ں���������ɫӳ��
}

static thread_local Array<NeighborRecord> neighs;

SummedDensityColorizer::SummedDensityColorizer(const RunSettings& settings, Palette palette)//�ܶȹ�����ɫ�����캯��
    : palette(std::move(palette)) {
    finder = Factory::getFinder(settings);//�����ܶȵ������finder
    kernel = Factory::getKernel<3>(settings);//����3ά�˺��������kernel
}

void SummedDensityColorizer::initialize(const Storage& storage, const RefEnum ref) {//�ܶȹ�����ɫ����ʼ��
    //��storage��ȡ��������λ�����ݵ�m��r
    m = makeArrayRef(storage.getValue<Float>(QuantityId::MASS), ref);
    r = makeArrayRef(storage.getValue<Vector>(QuantityId::POSITION), ref);

    finder->build(SEQUENTIAL, r);//��λ������r�����ܶȵ������finder
}

float SummedDensityColorizer::sum(const Size idx) const {//�ۼ��ܶȹ���ֵ�������ܶȺ�
    finder->findAll(idx, r[idx][H] * kernel.radius(), neighs);// ��������ֵidx�ڰ뾶��Χ���ҵ��������������ڵ�
    Float rho = 0._f;// ��ʼ���ܶȺ�Ϊ0
    for (const auto& n : neighs) {// �����������ڵ�
        rho += m[n.index] * kernel.value(r[idx] - r[n.index], r[idx][H]);// �ۼ�ÿ�����ڵ���ܶȹ���ֵ
    }
    return float(rho);
}

bool TemperatureColorizer::hasData(const Storage& storage) const {// ���洢�����Ƿ������¶���ɫ������
    if (!storage.has(QuantityId::DENSITY) || !storage.has(QuantityId::ENERGY)) {// // ����Ƿ����ܶȺ���������
        return false;
    }
    if (storage.getMaterialCnt() == 0) {// ����Ƿ��в�������
        return false;
    }
    for (Size matId = 0; matId < storage.getMaterialCnt(); ++matId) { // �������Ƿ���EosMaterial����
        const IMaterial& material = storage.getMaterial(matId);
        if (!dynamic_cast<const EosMaterial*>(&material)) {
            return false;
        }
    }
    return true;
}

void TemperatureColorizer::initialize(const Storage& storage, const RefEnum UNUSED(ref)) {//�¶���ɫ����ʼ��
    // ��ȡ�ܶȺ���������
    ArrayView<const Float> rho = storage.getValue<Float>(QuantityId::DENSITY);
    ArrayView<const Float> u = storage.getValue<Float>(QuantityId::ENERGY);
    Array<Float> T(rho.size());// ��ȡ�ܶȺ���������
    for (Size matId = 0; matId < storage.getMaterialCnt(); ++matId) {// ��ÿ������
        MaterialView view = storage.getMaterial(matId); // ��ȡ������ͼ
        EosMaterial& material = dynamic_cast<EosMaterial&>(view.material()); // ��ȡ���϶���
        for (Size i : view.sequence()) {// �Բ����е�ÿ����
            T[i] = material.getEos().getTemperature(rho[i], u[i]);// �����¶�
        }
    }
    values = makeArrayRef(std::move(T));//��Ҫ��ɫ��ֵ�����¶�����
}

void CorotatingVelocityColorizer::initialize(const Storage& storage, const RefEnum ref) {//��ת�ٶ���ɫ����ʼ��
    // ��ȡλ�á��ٶȺͲ���id
    r = makeArrayRef(storage.getValue<Vector>(QuantityId::POSITION), ref);
    v = makeArrayRef(storage.getDt<Vector>(QuantityId::POSITION), ref);
    matIds = makeArrayRef(storage.getValue<Size>(QuantityId::MATERIAL_ID), ref);

    ArrayView<const Float> m = storage.getValue<Float>(QuantityId::MASS);// ��ȡ����
    // ��λ�����ݺ��ٶ�����ת��Ϊ������ͼ
    ArrayView<const Vector> rv = r;
    ArrayView<const Vector> vv = v;
    data.resize(storage.getMaterialCnt());// ��ʼ��������������

    for (Size i = 0; i < data.size(); ++i) {// Ϊÿ�����ϼ������ĺͽ��ٶ�
        MaterialView mat = storage.getMaterial(i);// ��ȡ������ͼ
        // ����ò������ӵ�������Χ
        const Size from = *mat.sequence().begin();
        const Size to = *mat.sequence().end();
        const Size size = to - from;// ����ò��ϵ���������
        //  ʹ��getCenterOfMass��������ò��ϵ�����
        data[i].center = Post::getCenterOfMass(m.subset(from, size), rv.subset(from, size));
        //  ʹ��getAngularFrequency��������ò��ϵĽ��ٶ�
        data[i].omega =
            Post::getAngularFrequency(m.subset(from, size), rv.subset(from, size), vv.subset(from, size));
    }
}

bool DamageActivationColorizer::hasData(const Storage& storage) const {// ���洢�����Ƿ��������˼�����ɫ������
    // ����Ƿ������������:
    // 1. Ӧ��������ƫ΢��(Deviatoric stress)
    // 2. ѹǿ(Pressure)
    // 3. ��С��Ӧ��(Eps min)
    // 4. ����ֵ(Damage)
    return storage.has(QuantityId::DEVIATORIC_STRESS) && storage.has(QuantityId::PRESSURE) &&
           storage.has(QuantityId::EPS_MIN) && storage.has(QuantityId::DAMAGE);
}

void DamageActivationColorizer::initialize(const Storage& storage, const RefEnum UNUSED(ref)) {//���˼�����ɫ����ʼ��
    // ��ȡӦ������ƫ΢�֡�ѹǿ����С��Ӧ�䡢����ֵ
    ArrayView<const TracelessTensor> s = storage.getValue<TracelessTensor>(QuantityId::DEVIATORIC_STRESS);
    ArrayView<const Float> p = storage.getValue<Float>(QuantityId::PRESSURE);
    ArrayView<const Float> eps_min = storage.getValue<Float>(QuantityId::EPS_MIN);
    ArrayView<const Float> damage = storage.getValue<Float>(QuantityId::DAMAGE);

    ratio.resize(p.size());// ��ʼ����������
    /// \todo taken from ScalarGradyKippDamage, could be deduplicated ȡ��ScalarGradyKippDamage������ȥ��
    for (Size matId = 0; matId < storage.getMaterialCnt(); ++matId) { // ����ÿ������
        MaterialView mat = storage.getMaterial(matId);// ��ȡ������ͼ
        const Float young = mat->getParam<Float>(BodySettingsId::YOUNG_MODULUS); // ��ȡ����ģ��

        /// \todo parallelize ���л�
        for (Size i : mat.sequence()) {// ���������е�ÿ������
            const SymmetricTensor sigma = SymmetricTensor(s[i]) - p[i] * SymmetricTensor::identity();// �������ӵ�Ӧ������
            // ���Ӧ������������ֵ
            Float sig1, sig2, sig3;
            tie(sig1, sig2, sig3) = findEigenvalues(sigma);
            const Float sigMax = max(sig1, sig2, sig3);// �������ֵ
            const Float young_red = max((1._f - pow<3>(damage[i])) * young, 1.e-20_f);// �����С������ģ��
            const Float strain = sigMax / young_red;// �������Ӧ��
            ratio[i] = float(strain / eps_min[i]);// �������
        }
    }
}

BeautyColorizer::BeautyColorizer() {//������ɫ�����캯��
    palette = Palette({ { 0.f, Rgba(0.5f, 0.5f, 0.5) }, { 1.f, Rgba(1.f, 1.f, 0.6f) } },// ��ʼ����ɫ��
        Interval(u_0, u_yellow),
        PaletteScale::LOGARITHMIC);//�����߶�
    /// \todo maybe not fixed? �ܲ��ܲ��̶�
    // �ڵ�ɫ������������̶���ɫ��
    palette.addFixedPoint(u_glow, Rgba(0.5f, 0.5f, 0.5f));//����
    palette.addFixedPoint(u_red, Rgba(0.8f, 0.f, 0.f));//��ɫ
    f_glow = (log10(u_glow) - log10(u_0)) / (log10(u_yellow) - log10(u_0));//���㷢�������ӳ�����ӣ���������Ϊ[u_0, u_glow]��ͨ�������߶��µ����䳤�ȱ�������ӳ������
    //����ӳ��������Ϊ�˶Է�������������ӳ��, ʹ����������ɫ��������������֡�
}

BoundaryColorizer::BoundaryColorizer(const Detection detection, const Float threshold)//�߽���ɫ�����캯��
    : detection(detection) {
    if (detection == Detection::NEIGBOUR_THRESHOLD) {// ��������ֵ���:���н����ھӵ����ӱ���Ϊ�Ǳ߽����ӡ�������ֵ�ᵼ�±߽��������ࡣ
        neighbors.threshold = Size(threshold);//�趨��������ֵ
    } else {// ������ֵ���:�߽��ɽ������淨�ߵ����λ������ȷ���������ϸߣ��������е��������ھ����ӵ�������������ֵ�ᵼ�±߽����Ӽ��١�
        normals.threshold = threshold;//�趨������ֵ
    }
}

bool BoundaryColorizer::hasData(const Storage& storage) const {// ���洢�����Ƿ�����߽���ɫ������
    if (detection == Detection::NORMAL_BASED) {// �����ⷽʽΪ����
        return storage.has(QuantityId::SURFACE_NORMAL);// ����Ƿ���ڱ��淨������
    } else {// �����ⷽʽΪ������
        return storage.has(QuantityId::NEIGHBOR_CNT);// ����Ƿ��������������
    }
}

void BoundaryColorizer::initialize(const Storage& storage, const RefEnum ref) {//�߽���ɫ����ʼ��
    // ���ݼ�ⷽʽ��ȡ��Ӧ����
    if (detection == Detection::NORMAL_BASED) {// �����ⷽʽΪ����
        normals.values = makeArrayRef(storage.getValue<Vector>(QuantityId::SURFACE_NORMAL), ref);// ��ȡ���淨������
    } else {// �����ⷽʽΪ������
        neighbors.values = makeArrayRef(storage.getValue<Size>(QuantityId::NEIGHBOR_CNT), ref);// ��ȡ����������
    }
}

bool BoundaryColorizer::isInitialized() const {//�߽���ɫ����ʼ���ж�
    return (detection == Detection::NORMAL_BASED && !normals.values.empty()) ||//// ���ݼ�ⷽʽ����Ӧ����
           (detection == Detection::NEIGBOUR_THRESHOLD && !neighbors.values.empty());
}

Rgba BoundaryColorizer::evalColor(const Size idx) const {//�߽���ɫ���ķ�����ɫ
    if (isBoundary(idx)) {// ����isBoundary�����ж�����idx��Ӧ�������Ƿ��ڱ߽�
    //isBoundary�������ݹ��캯�����õı߽��ⷽʽ,ʹ����ֵ�ж�idx��Ӧ�������Ƿ��ڱ߽�
        return Rgba::red();//����Ǳ߽�����,��ɫΪ��ɫ
    } else {
        return Rgba::gray();//������Ǳ߽�����,��ɫΪ��ɫ
    }
}

bool BoundaryColorizer::isBoundary(const Size idx) const {//�߽��ж�
    switch (detection) {
    case Detection::NEIGBOUR_THRESHOLD:
        SPH_ASSERT(!neighbors.values.empty());//ȷ��������
        return neighbors.values[idx] < neighbors.threshold;//����������Ƿ������ֵ
    case Detection::NORMAL_BASED:
        SPH_ASSERT(!normals.values.empty());
        return getLength(normals.values[idx]) > normals.threshold;//��鷨�߳����Ƿ������ֵ
    default:
        NOT_IMPLEMENTED;
    }
}


/// \todo possibly move elsewhere ���ܻ�ᵽ�����ط�
static uint64_t getHash(const uint64_t value, const Size seed) {//����һ��64λ�޷��������Ĺ�ϣֵ
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

static Rgba getRandomizedColor(const Size idx, const Size seed = 0) {//����hash������������ԣ������������������������һ�������ɫֵ
    const uint64_t hash = getHash(idx, seed); //ʹ��getHash��������idx��seed����һ��64λ�޷��������Ĺ�ϣֵhash
    //��hash����ȡ����ɫ����ɫ����ɫ�����ֵ
    //r��ȡhash��16λ
    //g��ȡhash��16 - 32λ
    //b��ȡhash��32λ
    const uint8_t r = (hash & 0x00000000FFFF);
    const uint8_t g = (hash & 0x0000FFFF0000) >> 16;
    const uint8_t b = (hash & 0xFFFF00000000) >> 32;
    //��ȡ�������ֵӳ�䵽0 - 255�ķ�Χ��,ת����Ӧ�ĺ�����������ʹ�ú��������������һ��Rgba��ɫֵ�����ء�
    return Rgba(r / 255.f, g / 255.f, b / 255.f);
}

template <typename TDerived>//ģ�庯��
Rgba IdColorizerTemplate<TDerived>::evalColor(const Size idx) const {//IdColorizer����ɫ
    const Optional<Size> id = static_cast<const TDerived*>(this)->evalId(idx);//����evalId()������������ֵidx��ȡ��Ӧ��idֵ
    if (!id) {
        return Rgba::gray();//���idΪ��,ֱ�ӷ��ػ�ɫ
    }
    const Rgba color = getRandomizedColor(id.value(), seed);// ��id����,�����getRandomizedColor����id������seed���������ɫ
    return color;
}

template <typename TDerived>
Optional<Particle> IdColorizerTemplate<TDerived>::getParticle(const Size idx) const {//��ȡ����
    Particle particle(idx);//����һ�����Ӷ���
    const Optional<Size> id = static_cast<const TDerived*>(this)->evalId(idx);//��ȡid
    if (id) {
        particle.addValue(QuantityId::FLAG, id.value());//��id���ڣ����������ӵ�FLAG��������Ϊidxֵ
    }
    return particle;
}

//ģ����
template class IdColorizerTemplate<ParticleIdColorizer>;//����ID��ɫ��
template class IdColorizerTemplate<ComponentIdColorizer>;//���ID��ɫ��
template class IdColorizerTemplate<AggregateIdColorizer>;//��ȺID��ɫ��
template class IdColorizerTemplate<IndexColorizer>;//������ɫ��

void ParticleIdColorizer::initialize(const Storage& storage, const RefEnum ref) {//����ID��ɫ����ʼ��
    if (storage.has(QuantityId::PERSISTENT_INDEX)) {//���洢�д��ڳ�������
        persistentIdxs = makeArrayRef(storage.getValue<Size>(QuantityId::PERSISTENT_INDEX), ref);//�Ӵ洢�л�ȡ��������,���洢������
    }
}

Optional<Particle> ParticleIdColorizer::getParticle(const Size idx) const {//����ID��ɫ����ȡ����
    Particle particle(idx);
    particle.addValue(QuantityId::FLAG, idx);//�������ӵ�FLAG��������Ϊidxֵ
    if (!persistentIdxs.empty() && idx < persistentIdxs.size()) {//���־��������鲻Ϊ����������Ч
        particle.addValue(QuantityId::PERSISTENT_INDEX, persistentIdxs[idx]);//ȡ����Ӧ�־�����,��ӵ����ӵĶ�Ӧ������
    }
    return particle;
}


ComponentIdColorizer::ComponentIdColorizer(const GuiSettings& gui,//// ���ID��ɫ�����캯��
    const Flags<Post::ComponentFlag> connectivity,
    const Optional<Size> highlightIdx)
    : IdColorizerTemplate<ComponentIdColorizer>(gui)
    , connectivity(connectivity)//�洢�����������Ϣ
    , highlightIdx(highlightIdx) {}//�洢Ҫ�������������

void ComponentIdColorizer::setHighlightIdx(const Optional<Size> newHighlightIdx) {// ���ø����������
    if (newHighlightIdx) {
        highlightIdx = min(newHighlightIdx.value(), components.size() - 1);//�����������Ч,ȡ���������Ŀ����Сֵ��Ϊ��������
    } else {
        highlightIdx = NOTHING;//���������Ϊ��,�����ø�������Ϊ��Чֵ
    }
}

Rgba ComponentIdColorizer::evalColor(const Size idx) const {//���ID��ɫ����ɫ
    if (highlightIdx) {// ������ڸ�������
        if (highlightIdx.value() == components[idx]) {//����ǰ���Ϊ�������
            return Rgba(1.f, 0.65f, 0.f);//������ɫ
        } else {
            return Rgba::gray(0.3f);//������Ϊ������ɫ
        }
    } else { // û�и������� 
        const Size id = this->evalId(idx).value();// ��ȡ���id
        return getRandomizedColor(compIdxs[id], seed);// ����idѡ�������ɫ
    }
}

Optional<Particle> ComponentIdColorizer::getParticle(const Size idx) const {// ���ID��ɫ����ȡ����
    Particle particle(idx);// �������Ӷ���
    const Optional<Size> id = this->evalId(idx); // ��ȡ���id
    particle.addValue(QuantityId::FLAG, id.value());//�������ӵ�FLAG����
    // �ҵ�����ͬһ�����������������
    Array<Size> indices;
    for (Size i = 0; i < r.size(); ++i) {
        if (components[i] == id.value()) {
            indices.push(i);
        }
    }
    if (indices.size() > 1) { // ���������������,��������Ľ��ٶ�
        const Vector omega = Post::getAngularFrequency(m, r, v, indices);
        particle.addValue(QuantityId::ANGULAR_FREQUENCY, getLength(omega));//�������ӵ�ANGULAR_FREQUENCY����
    }
    return particle;
}

bool ComponentIdColorizer::hasData(const Storage& storage) const {// ���洢�����Ƿ��������ID��ɫ������
    return hasVelocity(storage);// ���洢���Ƿ�����ٶ�����
}

void ComponentIdColorizer::initialize(const Storage& storage, const RefEnum ref) {//���ID��ɫ����ʼ��
    const Array<Vector>& current = storage.getValue<Vector>(QuantityId::POSITION); //��ȡ��ǰλ������
    if (current == cached.r) {//�ͻ���λ�����ݱȽϣ������ͬ,������ʼ��;�����ͬ,����г�ʼ��
        // optimization, very poorly done �Ż������ĺܲ�
        return;
    }

    m = makeArrayRef(storage.getValue<Float>(QuantityId::MASS), ref);//��ȡ��������
    r = makeArrayRef(storage.getValue<Vector>(QuantityId::POSITION), ref);//��ȡλ������
    v = makeArrayRef(storage.getDt<Vector>(QuantityId::POSITION), ref);//��ȡ�ٶ�����

    cached.r = current.clone();//����ǰλ������current���ƴ洢������cached.r��

    const Size numComponents = Post::findComponents(storage, 2.5_f, connectivity, components);// ʹ��findComponents���������������

    // sort by the smallest index in each component ��ÿ���������С����������
    //��ʼ��compIdxs����ΪnumComponents��С, ���Ϊ-1
    compIdxs.resize(numComponents);
    compIdxs.fill(Size(-1));
    for (Size i = 0; i < components.size(); ++i) {//����components���飬��¼�������С��������
        const Size ci = components[i];
        compIdxs[ci] = min(compIdxs[ci], i);
    }
}

String ComponentIdColorizer::name() const {//����connectivity�����÷��ز�ͬ������,�Է�ӳ����Ļ��ַ�ʽ
    if (connectivity.has(Post::ComponentFlag::ESCAPE_VELOCITY)) {//������Ӳ������Ƿ�������ESCAPE_VELOCITY���������ٶȻ��֣�����
        return "Bound component ID";//����"Bound component ID"
    } else if (connectivity.has(Post::ComponentFlag::SEPARATE_BY_FLAG)) {//������Ӳ������Ƿ�������SEPARATE_BY_FLAG����FLAG���֣�����
        return "Component ID (flag)";//����"Component ID (flag)"
    } else {//�����δ����
        return "Component ID";//����Ĭ��"Component ID"
    }
}

void MaterialColorizer::initialize(const Storage& storage, const RefEnum ref) {//������ɫ����ʼ��
    IndexColorizer::initialize(storage, ref);// ���ø����initialize��������ʼ��

    const Size matCnt = storage.getMaterialCnt();// ��ȡ��������
    //��ʼ���洢EOS��Rheology���͵�����
    //EOS��Equation of state for this material ״̬����
    //Rheology��Model of stress reducing used within the rheological model ����ģ����ʹ�õ�Ӧ����Сģ��
    eosNames.resize(matCnt);
    rheoNames.resize(matCnt);
    for (Size matId = 0; matId < matCnt; ++matId) {//����ÿ������
        const IMaterial& mat = storage.getMaterial(matId);//��ȡ������Ϣ
        //��ȡEOS��Rheology���ͣ��������Ӧ����
        const EosEnum eos = mat.getParam<EosEnum>(BodySettingsId::EOS);
        const YieldingEnum yield = mat.getParam<YieldingEnum>(BodySettingsId::RHEOLOGY_YIELDING);
        eosNames[matId] = EnumMap::toString(eos);
        rheoNames[matId] = EnumMap::toString(yield);
    }
}

Optional<Particle> MaterialColorizer::getParticle(const Size idx) const {//������ɫ����ȡ����
    Particle particle(idx);
    const Optional<Size> id = IndexColorizer::evalId(idx);// // ��ȡ����id
    if (id) {//id����
        particle.addValue(QuantityId::MATERIAL_ID, id.value());//�������ӵ�MATERIAL_ID����
        //�������ӵ�EOS��RHEOLOGY_YIELDING����
        particle.addParameter(BodySettingsId::EOS, eosNames[id.value()]);
        particle.addParameter(BodySettingsId::RHEOLOGY_YIELDING, rheoNames[id.value()]);
    }
    return particle;
}

void TimeStepColorizer::initialize(const Storage& storage, const RefEnum ref) {//ʱ�䲽����ɫ����ʼ��
    TypedColorizer<Float>::initialize(storage, ref);//���ø����initialize��������ʼ��
    critIds = makeArrayRef(storage.getValue<Size>(QuantityId::TIME_STEP_CRITERION), ref);//��ȡʱ�䲽���ж���׼
}

Optional<Particle> TimeStepColorizer::getParticle(const Size idx) const {//ʱ�䲽����ɫ����ȡ����
    Particle particle(idx);
    particle.addValue(QuantityId::TIME_STEP, values[idx]);// �������ӵ�ʱ�䲽������
    particle.addValue(QuantityId::TIME_STEP_CRITERION, toString(CriterionId(critIds[idx])));// �������ӵ�ʱ�䲽���ж���׼����
    return particle;
}

NAMESPACE_SPH_END
