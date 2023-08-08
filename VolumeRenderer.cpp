#include "gui/renderers/VolumeRenderer.h"
#include "gui/Factory.h"
#include "gui/objects/Camera.h"
#include "gui/objects/Colorizer.h"
#include "gui/renderers/FrameBuffer.h"
#include "objects/finders/KdTree.h"
#include "objects/utility/OutputIterators.h"

NAMESPACE_SPH_BEGIN
//�����Ⱦ�������Ⱦ��һ��������ʾ3D��ɢ�������ݼ���ͨ����3D����������2DͶӰ�ļ����������֮�����ǽ��������أ�����Ļ���Ԫ�أ���������ȾΪ�ɼ���2Dͼ��
//����Ӿ�Ч�������������ơ����棬�Ǻ����ü�������ֵġ���ʱ�����ǿ���ʹ�������Ⱦ��������ͨ���������������еĴ�����ʽ���໥���ã����ɱ�������Ч����
VolumeRenderer::VolumeRenderer(SharedPtr<IScheduler> scheduler, const GuiSettings& settings)//���������Ⱦ����
    : IRaytracer(scheduler, settings) {}//�̳��Թ�׷��

VolumeRenderer::~VolumeRenderer() = default;//��������

const float MAX_DISTENTION = 50;//������ͳ���
const float MIN_NEIGHS = 8;//��С�ھӳ���

void VolumeRenderer::initialize(const Storage& storage, const IColorizer& colorizer, const ICamera& camera) {//��ʼ��������Storage��IColorizer��ICamera��Ϊ����
    cached.r = storage.getValue<Vector>(QuantityId::POSITION).clone();//�Ӵ洢�л�ȡλ������,����¡��cached.r��
    //���ﻺ��洢���ǵ��λ����Ϣ��
    this->setColorizer(colorizer);//������ɫ�������ڽ�λ������ӳ��ΪRGB��ɫ��

    cached.distention.resize(cached.r.size());//�洢���ӵ�����ϵ����

    KdTree<KdNode> tree;
    //�����ǳ������д��������ʱ��ÿ�����߶�Ҫ��ÿ�������󽻵������Ƿǳ�����ģ�������Ҫ����������̽��м��١�
    //���ٷ���ͨ����Ϊ�ռ仮�֣����廮�֡��ռ仮���а˲�����KDTree,BSPTree,���廮����BVH��SAH�ȡ�
    tree.build(*scheduler, cached.r, FinderFlag::SKIP_RANK);//����һ��kd�������������������־,���Լ��ٹ������̣���
    //kd����k-dimensional���ļ��,����KdTreeʱ,��ͨ���ݹ�ָ�ķ�ʽ,���㼯����Ϊ����ڵ�,ÿ���ڵ������ٻ����ӽڵ�,ֱ���ﵽ�����ȡ�
    //���������һ��������νڵ�����ṹ, ͨ�����ֿռ仮��,�������Լ����ҵ��ڽ���Ĳ�ѯ��
    //���ں������ھӵ���ң����Լ��ٲ���λ����Χ�ĵ�, �Ӷ��Ż���ȾЧ�ʡ�
    Array<BvhSphere> spheres(cached.r.size());//���ڴ洢ÿ����İ�Χ��Χ��
    spheres.reserve(cached.r.size());
    ThreadLocal<Array<NeighborRecord>> neighs(*scheduler);//ʹ�� neighs �߳̽��б��ش洢,Ϊÿ���߳�ά�������洢��
    parallelFor(*scheduler, neighs, 0, cached.r.size(), [&](const Size i, Array<NeighborRecord>& local) {//����ÿ����
        const float initialRadius = cached.r[i][H];//�����ʼ�뾶
        float radius = initialRadius;
        while (radius < MAX_DISTENTION * initialRadius) {//�� initialRadius ��ʼ����뾶�����ھӵ㡣
            tree.findAll(cached.r[i], radius, local); //ʹ�� KdTree �� findAll �������ҵ�ǰ�뾶�µ��ھӵ㡣
            if (local.size() >= MIN_NEIGHS) {//����ھӵ��� �� MIN_NEIGHS,��ȷ����ǰ�뾶Ϊ���հ뾶��
                break;
            } else {
                radius *= 1.5f;//����,���뾶����1.5��,����������
            }
        }

        BvhSphere s(cached.r[i], radius);//������Χ�����ݣ��뾶Ϊȷ�������հ뾶��
        s.userData = i;
        spheres[i] = s;// ����Χ��д�� spheres ���顣

        cached.distention[i] = min(radius / initialRadius, MAX_DISTENTION);//�����������ϵ����ȡ���հ뾶���ʼ�뾶�ı�ֵ��������ͳ����е���Сֵ��
    });

    //���������Ͳ�����Ϣ�����ÿ����ĺ���ο��뾶
    ArrayView<const Float> m = storage.getValue<Float>(QuantityId::MASS);// ��ȡ����������m��
    cached.referenceRadii.resize(cached.r.size());//�����ο��뾶����Ĵ�С����������ͬ��
    if (storage.getMaterialCnt() > 0) {
        for (Size matId = 0; matId < storage.getMaterialCnt(); ++matId) {//����ÿ������
            MaterialView mat = storage.getMaterial(matId);
            const Float rho = mat->getParams().has(BodySettingsId::DENSITY)//��ȡ�����ܶ�rho��
                                  ? mat->getParam<Float>(BodySettingsId::DENSITY)
                                  : 1000._f;
            for (Size i : mat.sequence()) {//�����ڸò��ʵĵ�i:
                const Float volume = m[i] / rho;//��������� = ���� / �ܶ�
                cached.referenceRadii[i] = root<3>(3._f * volume / (4._f * PI));//ͨ���������ο��뾶
            }
        }
    } else {//���û�в�����Ϣ
        // guess the dentity
        const Float rho = 1000._f;//ʹ��Ĭ���ܶ�1000������ο��뾶
        for (Size i = 0; i < m.size(); ++i) {
            const Float volume = m[i] / rho;
            cached.referenceRadii[i] = root<3>(3._f * volume / (4._f * PI));
        }
    }

    cached.textures.clear();//�����������
    for (Size i = 0; i < storage.getAttractorCnt(); ++i) {//��������������attractor
        const Attractor& a = storage.getAttractors()[i];//��storage��ȡ����attractor
        const bool visible = a.settings.getOr(AttractorSettingsId::VISIBLE, true);//��ȡattractor�Ŀɼ��Բ���
        const Float albedo = a.settings.getOr(AttractorSettingsId::ALBEDO, 1._f);//��ȡattractor�ķ����ʲ���
        cached.attractors.push(AttractorData{ a.mass, a.position, a.radius, visible, albedo });//�������Ӳ�������cached.attractors���顣

        String texturePath = a.settings.getOr<String>(AttractorSettingsId::VISUALIZATION_TEXTURE, "");//��ȡattractor������·��
        if (!texturePath.empty()) {//���������
            if (cached.textureCache.contains(texturePath)) {//����������Ѵ���,ֱ��push��cached.textures����
                cached.textures.push(cached.textureCache[texturePath]);
            } else {//�����������,���뻺�沢push��
                SharedPtr<Texture> texture =
                    makeShared<Texture>(Path(texturePath), TextureFiltering::BILINEAR);
                cached.textureCache.insert(texturePath, texture);
                cached.textures.push(texture);
            }
        } else {// ���������,pushһ��nullptr
            cached.textures.push(nullptr);
        }

        BvhSphere sphere(a.position, a.radius);//����attractor�İ�Χ��,����spheres���顣
        sphere.userData = cached.r.size() + i;
        spheres.push(sphere);
    }

    bvh.build(std::move(spheres));//.ʹ�ð�Χ�����ݹ������ٽṹbvh��
    //��ΰ�Χ��BVH��Bounding Volume Hierachies����BV��Bounding Volume���ǰ���һ������Ŀռ��壬�������������ļ���������״Ҫ�򵥵Ķ࣬���Զ���������ײ����ٶȱ�ֱ�������屾�����ཻ���졣
    //BVH��һ��������BVΪ�������л��ֵĽṹ����������һ���������ɸ��ڵ㡢�ڲ��ڵ��Ҷ�ӽڵ���ɡ�����Ҷ�ӽڵ������壬ÿ����Ҷ�ӽڵ㶼�а�Χ�壬���ڵ���԰��ӽڵ��Χ������
    
    //���㳡��������
    cached.maxDistance = 0;
    for (const Attractor& a : storage.getAttractors()) {//��������attractor
        const Float dist = getLength(a.position - camera.getPosition());//��ÿ��attractor,������������ľ���dist��
        cached.maxDistance = max(cached.maxDistance, 2 * dist);//��dist��2���뵱ǰ��cached.maxDistance�Ա�,�����ϴ�ֵ��Ϊ�µ�cached.maxDistance��
    }

    for (ThreadData& data : threadData) {//��ʼ���߳�����threadData
        data.data = RayData{};
    }

    shouldContinue = true;//������Ⱦ״̬shouldContinueΪtrue��ȷ����Ⱦ���Կ�ʼ
}

bool VolumeRenderer::isInitialized() const {//��黺�������Ƿ�Ϊ�գ��Դ��ж��Ƿ�����˳�ʼ����
    return !cached.r.empty();
}

void VolumeRenderer::setColorizer(const IColorizer& colorizer) {//������ɫ��
    cached.colors.resize(cached.r.size());//������ɫ�����С
    for (Size i = 0; i < cached.r.size(); ++i) {//����evalColor������ɫ
        cached.colors[i] = colorizer.evalColor(i);
    }
}

Rgba VolumeRenderer::shade(const RenderParams& params, const CameraRay& cameraRay, ThreadData& data) const {//��Ӱ��ɫ
    const Vector primaryDir = getNormalized(cameraRay.target - cameraRay.origin);//������������ߵ�������primaryDir:��originָ��target
    const Ray primaryRay(cameraRay.origin, primaryDir);//ʹ�������򹹽�������primaryRay

    RayData& rayData(data.data);//�߳����ݣ���������Ⱦ����
    LensingEffect::Segments& segments = rayData.segments;//���߶�����
    Array<CurvedRayIntersectionInfo>& intersections = rayData.intersections;//��������

    /// ����LensingEffect��,����͸��ЧӦ��ȡ���еĽ��㡣
    LensingEffect lensing(cached.attractors,
        params.relativity.lensingMagnitude,
        0.1_f * cached.maxDistance,
        cached.maxDistance,
        params.volume.absorption > 0.f);
    //����:attractors��lensing���������������롢����ϵ��������lensing����
    //���: ���߶�segments�ͽ���intersections
    const Ray lastRay = lensing.getAllIntersections(bvh, primaryRay, segments, intersections);
    //����LensingEffect��getAllIntersections����,����������primaryRay��������͸��ЧӦ�����н�������߶Ρ��洢��segments��intersections�С�
    //�õ���������lastRay
    Rgba result = this->getEnviroColor(CameraRay{ lastRay.origin(), lastRay.origin() + lastRay.direction() });//����lastRay���㻷������ɫ

    for (const CurvedRayIntersectionInfo& is : reverse(intersections)) {// ����������й��߽���
        const BvhSphere* s = static_cast<const BvhSphere*>(is.object);// ���������ת��Ϊ����ָ��
        const Size i = s->userData;// ��ȡ�����userData��Ϊ��������
        const Vector hit = is.segment->origin() + is.segment->direction() * is.t; // ������ײ������ = ����ԭ�� + ���� * �������
        const Vector center = s->getCenter();// ��ȡ��������
        const Vector toCenter = getNormalized(center - hit);// ������ײ��ָ�����ķ����������������һ����Ϊ��λ����
        
        //�����������attractor�����
        //��������һ����������ӣ�������������򣬶�һ�����Ӳ���Ӱ�����á�
        if (i >= cached.r.size()) {//// �����������������Ŀ,��Ϊattractor
            // attractor, a solid object -> erase the emission accumulated so far if visible  ��������һ����壬��Ҫͨ����ǰ��ɫ�����Ŀǰ�ۻ��Ŀɼ��ķ�����
            const Size idx = i - cached.r.size();// ���������� = �������� - ������Ŀ
            if (cached.attractors[idx].visible) {//����ɼ�
                result = this->getAttractorColor(params, idx, hit);//����getAttractorColor����attractor����ײ��hit����ɫ
            }
            continue;
        }
        // ����������,���б�׼��ɫ����
        const float cosPhi = abs(dot(toCenter, is.segment->direction()));// ������ײ�㵽������������߷��������ֵ
        const float distention = cached.distention[i]; //��ȡ��i�����ӵ�����ϵ��
        // smoothing length should not have effect on the total emission ƽ�����Ȳ�Ӧ���ܷ���������Ӱ��
        const float radiiFactor = cached.referenceRadii[i] / cached.r[i][H];//�������Ӱ뾶����radiiFactor,�òο��뾶���Ե�ǰ�뾶,�������⻬���ȵ�Ӱ�졣
        const float secant = 2._f * getLength(center - hit) * cosPhi * radiiFactor;//���㽻�㵽���ĵĸ��߳���secant,���ں������㷢��ǿ��
        // make dilated particles absorb more ʹ���͵��������ո������
        result = result * exp(-params.volume.absorption * secant * pow<3>(cosPhi/distention));// ʹ�� secant �� cosPhi ����������,���������ͳ̶ȵ��������ǿ
        // 3th power of cosPhi to give more weight to the sphere center,
        // divide by distention^3; distention should not affect the total emission  ����cospi��3�η���Ϊ��ʹ���ĵ�Ӱ��Ȩ�ظ��󣻳�������ϵ�������η���Ϊ��ʹ���Ͳ�Ӱ���ܷ�������
        const float magnitude = params.volume.emission * pow<3>(cosPhi / distention) * secant;// ���㷢����������Ȩ�ظ�,��������Ӱ��
        result += cached.colors[i] * magnitude;//������������������ɫ�õ������ӵķ���Ч������������ӵ�����С�
        result.a() += magnitude;//����������ӵ�alphaͨ����������ͨ��(Alpha Channel)��ָһ��ͼƬ��͸���Ͱ�͸����
    }
    result.a() = min(result.a(), 1.f);// ������alpha������1����
    return result;// ���ؽ��
}

Rgba VolumeRenderer::getAttractorColor(const RenderParams& params,//��������ɫ
    const Size index,
    const Vector& hit) const {
    const AttractorData& a = cached.attractors[index];//// ��ȡ��ǰattractor����
    Rgba diffuse = Rgba::gray(a.albedo);// ��ʼ����������ɫΪ��ɫ
    const SharedPtr<Texture>& texture = cached.textures[index];//// ��ȡattractor����
    if (texture) {//���������
        const Vector r0 = hit - a.position;//������ײ����attractor����ϵ�µ���������
        SphericalCoords spherical = cartensianToSpherical(r0);//����������ת��Ϊ������
        Vector uvw = Vector(0.5_f - spherical.phi / (2._f * PI), spherical.theta / PI, 0._f);// ��������ת��Ϊuvw��������
        diffuse = texture->eval(uvw) * a.albedo;//����������albedo������ϵ�������Ϊ������������ɫ
    }

    const Vector n = getNormalized(a.position - hit);//������ײ��ָ��attractor���򣬾�����һ����Ϊ��λ����
    const Float cosPhi = dot(n, params.lighting.dirToSun);//����n��������շ���нǵ�����ֵcosPhi
    if (cosPhi <= 0._f) {// �������<=0,������n���������ֱ 
        // not illuminated -> just ambient light  ���ᱻ�������䵽,�����ػ�����
        return diffuse * params.lighting.ambientLight;//���ػ�������ɫ
    }
        return diffuse * float(PI * cosPhi * params.lighting.sunLight + params.lighting.ambientLight);//���򷵻�: ��������ɫ*(�������� + ̫������ * ����)
   }

NAMESPACE_SPH_END
