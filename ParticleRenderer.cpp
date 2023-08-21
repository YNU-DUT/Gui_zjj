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
//������Ⱦ��������Ⱦ����Ч��
static void drawVector(IRenderContext& context,//����һ������ͷ�߶α�ʾ����
    const ICamera& camera,
    const Vector& r,
    const Vector& v,
    const Float length) {
    if (getSqrLength(v) == 0._f) {// �������v����Ϊ0,������
        return;
    }
    const Optional<ProjectedPoint> p1 = camera.project(r); // ͶӰ�������
    const Optional<ProjectedPoint> p2 = camera.project(r + v);// ͶӰ�յ�(���+��������)
    if (!p1 || !p2) {//ֻҪp1��p2����һ����Ч����ֹͣ���ơ�
        return;
    }

    Coords dir = p2->coords - p1->coords;//����p1ָ��p2���������dir,
    const Float l = getLength(dir);//��ȡ������������
    if (l == 0._f) {//�������Ϊ0���򲻻���
        return;
    }
    dir *= float(length / l);// ��dir���ŵ�ָ������ length
    const Coords c1 = p1->coords;// �������
    const Coords c2 = p1->coords + dir;// �յ�����(���+��������)

    context.setColor(Rgba(1.f, 0.65f, 0.f), ColorFlag::LINE);// �����߶���ɫ
    context.setThickness(2.f);// �����߶δ�ϸ
    context.drawLine(c1, c2);// �����߶�

    // make an arrow ����ͷ
    AffineMatrix2 rot = AffineMatrix2::rotate(160._f * DEG_TO_RAD);// ����һ��˳ʱ����ת160�ȵı任����
    PlotPoint dp(dir.x, dir.y);// ����������ת��Ϊ��,�Խ��о���任
    PlotPoint a1 = rot.transformPoint(dp) * 0.1f;// ��dp��Ӧ����ת����,�õ���ͷһ���ƫ����
    PlotPoint a2 = rot.transpose().transformPoint(dp) * 0.1f;// ��dp��Ӧ��ת����ת����,�õ���ͷ��һ���ƫ���� 

    context.drawLine(c2, c2 + Coords(float(a1.x), float(a1.y)));// ���յ������ָ��a1���߶�,�γɼ�ͷ��һ��
    context.drawLine(c2, c2 + Coords(float(a2.x), float(a2.y)));// ���յ������ָ��a1���߶�,�γɼ�ͷ��һ��
}

static void drawGrid(IRenderContext& context, const ICamera& camera, const float grid) {//����һ���������������񱳾�
    // find (any) direction in the camera plane �����ƽ�����ҵ�һ������
    const Optional<CameraRay> originRay = camera.unproject(Coords(0, 0));// ����Ļ���ķ���һ������
    const Vector dir = getNormalized(originRay->target - originRay->origin);// �������߷�������dir,����һ��
    // ������dir��ֱ������perpDir
    Vector perpDir;
    if (dir == Vector(0._f, 0._f, 1._f)) {// ���dir��ֱ����
        perpDir = Vector(1._f, 0._f, 0._f);//��perpDirȡ(1,0,0)
    } else {
        perpDir = getNormalized(cross(dir, Vector(0._f, 0._f, 1._f)));// ����ȡdir��(0,0,1)�Ĳ����ΪperpDir
    }

    // find how much is projected grid distance ���ͶӰ���������
    const Coords shifted = camera.project(originRay->origin + grid * perpDir)->coords;//����ԭ���ƶ�һ�����񳤶�����Ļ�ϵ�ͶӰƫ����
    const float dx = getLength(shifted);//ȡdx��ֵΪshifted����Ϊ��ֱ�����ߵ�λ���
    const float dy = dx;// dyȡ��dx��ͬ��С����Ϊˮƽ�����ߵ�λ���
    const Coords origin = camera.project(Vector(0._f))->coords;// ԭ������ 

    context.setColor(Rgba(0.16f), ColorFlag::LINE);// ���û�����ɫ
    const Pixel size = context.size();// ���û�����ɫ
    for (float x = origin.x; x < size.x; x += dx) {// ��x�����ϣ���origin.x��ʼ,����dx,���ƴ�ֱ������
        context.drawLine(Coords(x, 0), Coords(x, size.y));
    }
    for (float x = origin.x - dx; x >= 0; x -= dx) {//��ͬ�������ڷ�������ѭ������
        context.drawLine(Coords(x, 0), Coords(x, size.y));
    }
    for (float y = origin.y; y < size.y; y += dy) {// ��y�����ϣ���origin.y��ʼ,����dy,����ˮƽ������
        context.drawLine(Coords(0, y), Coords(size.x, y));
    }
    for (float y = origin.y - dy; y >= 0; y -= dy) {//��ͬ�������ڷ�������ѭ������
        context.drawLine(Coords(0, y), Coords(size.x, y));
    }
}

ParticleRenderer::ParticleRenderer(const GuiSettings& settings) {//ParticleRenderer�Ĺ��캯��
    grid = float(settings.get<Float>(GuiSettingsId::VIEW_GRID_SIZE)); // ��GuiSettings������,��ȡVIEW_GRID_SIZE����ָ����������ֵ
  // ��ֵ�����ں����������񱳾�ʱ���������ߵļ��
    shouldContinue = true;// ��ʼ��shouldContinue��־Ϊtrue
  // �ñ�־���ڿ�����Ⱦѭ���ļ���ִ��
}

static bool isCutOff(const Vector& r, const Optional<float> cutoff, const Vector direction) {//�ж�һ�������Ƿ��ڽضϰ뾶��
    return cutoff && abs(dot(direction, r)) > cutoff.value();//���ָ����cutoff,��r��direction�����ϵ�ͶӰ���ȵľ���ֵ����cutoff,�򷵻�true,��ʾ�ڽضϰ뾶�⡣
    //���򷵻�false,��ʾ�ڽضϰ뾶�ڡ�
}

constexpr Size GHOST_INDEX = Size(-1);//������������
constexpr Size ATTRACTOR_INDEX = Size(-2);//����������

void ParticleRenderer::initialize(const Storage& storage,//������Ⱦ����ʼ��
    const IColorizer& colorizer,
    const ICamera& camera) {
    MEASURE_SCOPE("ParticleRenderer::initialize");// ��¼��������ʱ��
    cached.idxs.clear();// ���������������
    cached.positions.clear(); // �������λ�û��� 
    cached.colors.clear(); // ���������ɫ����
    cached.vectors.clear();// ���������ɫ����������

    //������������,ɸѡ�ڽضϷ�Χ�ڵ�,��ǰ��������Ļ�ռ���Ϣ,Ϊ������Ⱦ��׼����
    const Optional<float> cutoff = camera.getCutoff();// ��ȡ����Ľضϰ뾶cutoff
    const Vector direction = camera.getFrame().row(2);//��ȡ������ӽǷ���direction
    bool hasVectorData = bool(colorizer.evalVector(0));// �ж���ɫ���Ƿ�����������
    ArrayView<const Vector> r = storage.getValue<Vector>(QuantityId::POSITION);//��ȡ�������ӵ�λ������r
    for (Size i = 0; i < r.size(); ++i) {//����ÿ������r[i]
        const Optional<ProjectedPoint> p = camera.project(r[i]);// ������λ��ͶӰ����Ļ
        if (p && !isCutOff(r[i], cutoff, direction)) {// ����ڽضϰ뾶�������
            cached.idxs.push(i);// ������������
            cached.positions.push(r[i]);//��������λ��

            const Rgba color = colorizer.evalColor(i);//����������ɫ
            cached.colors.push(color);//����������ɫ

            if (hasVectorData) {//�����ɫ������������
                Optional<Vector> v = colorizer.evalVector(i);//��ȡ�������ӵ���ɫ����������
                SPH_ASSERT(v);//�����������ݲ�Ϊ��
                cached.vectors.push(v.value());//������ɫ����������
            }
        }
    }

    //������������
    SharedPtr<IStorageUserData> data = storage.getUserData();// �Ӵ洢�л�ȡ�Զ����û����� 
    if (RawPtr<GhostParticlesData> ghosts = dynamicCast<GhostParticlesData>(data.get())) {// �ж��Ƿ���GhostParticlesData���ͣ��Ӷ��ж��ǲ�����������
        for (Size i = 0; i < ghosts->size(); ++i) {// ����ÿ����������
            const Vector pos = ghosts->getGhost(i).position;// ��ȡ��������λ����Ϣ
            if (!isCutOff(pos, cutoff, direction)) {// ���������ضϰ뾶��
                cached.idxs.push(GHOST_INDEX);// ��GHOST_INDEX��Ϊ������������
                cached.positions.push(pos);//�����������ӵ�λ����Ϣ
                cached.colors.push(Rgba::transparent());//������ɫ����Ϊ͸��

                if (hasVectorData) {//�����ɫ������������
                    cached.vectors.push(Vector(0._f));//��������������Ϊ0
                }
            }
        }
    }

    //����������
    for (const Attractor& a : storage.getAttractors()) {//��ȡ�洢�е����������Ӳ�����
        if (!isCutOff(a.position, cutoff, direction)) {// ���������ضϰ뾶��
            cached.idxs.push(ATTRACTOR_INDEX);// ��ATTRACTOR_INDEX��Ϊ������������
            cached.positions.push(setH(a.position, a.radius));//���������ӵ�λ����Ϣ
            cached.colors.push(Rgba::white());//������ɫ����Ϊ��ɫ

            if (hasVectorData) {//�����ɫ������������
                cached.vectors.push(Vector(0._f));//��������������Ϊ0
            }
        }
    }

    // sort in z-order ��Z��˳�������������л������ݱ���ͬ������������Ⱦ������
    Order order(cached.positions.size());//����һ��Order����,���СΪ������λ����Ϣ�Ĵ�С
    order.shuffle([this, &direction](Size i, Size j) {// ʹ��shuffle�����Զ�Ԫν�ʱ������� 
        // �Ƚ����ӽǷ����ϵ�ͶӰ��С
        const Vector r1 = cached.positions[i];
        const Vector r2 = cached.positions[j];
        return dot(direction, r1) > dot(direction, r2);// �Ƚ������������ӽǷ����ϵ�ͶӰ��С
    });
    /// \todo could be changed to AOS to sort only once ���Ըĳ�AOS�㷨��ֻ����һ��
    // ��λ�á���������ɫӦ��������
    cached.positions = order.apply(cached.positions);
    cached.idxs = order.apply(cached.idxs);
    cached.colors = order.apply(cached.colors);

    cached.cameraDir = direction;// ��λ�á���������ɫӦ��������

    if (hasVectorData) {// �����ɫ������������
        cached.vectors = order.apply(cached.vectors);// ����������ҲӦ������
    } else {
        cached.vectors.clear();//���������������
    }
}

bool ParticleRenderer::isInitialized() const {//�ж��Ƿ��ʼ��
    return !cached.positions.empty();//ͨ���ж�cached.positions�Ƿ�Ϊ��,���ƶ�cached�е������Ƿ��ѱ���ʼ��
}

void ParticleRenderer::setColorizer(const IColorizer& colorizer) {//������ɫ��
    for (Size i = 0; i < cached.idxs.size(); ++i) {// ���������е�������������
        if (cached.idxs[i] == GHOST_INDEX ||
            cached.idxs[i] == ATTRACTOR_INDEX) { // ������������ӻ�������,������
            continue; // ghost or attractor
        }
        cached.colors[i] = colorizer.evalColor(cached.idxs[i]); // ������ͨ���ӵ���ɫ
    }
}

static AutoPtr<IRenderContext> getContext(const RenderParams& params, Bitmap<Rgba>& bitmap) {//��ȡ��Ⱦ�����ģ���Ⱦ��ǰ��׼���ӿ���
    if (params.particles.doAntialiasing) {// ������ÿ���ݣ�FSAA��������������ͼ�α�Ե�ľ�ݣ�ʹ���濴������Ϊƽ��
        if (params.particles.smoothed) {//�������ƽ����
            CubicSpline<2> kernel;// ʹ��CubicSpline���κ˺���
            return makeAuto<SmoothedRenderContext>(bitmap, kernel);// ����ƽ�����˵���Ⱦ������ 
        } else {
            return makeAuto<AntiAliasedRenderContext>(bitmap);// ���򴴽������������Ⱦ������
        }
    } else {//����ҪFSAA
        if (params.background.a() == 1.f) {// ���������͸��
            return makeAuto<PreviewRenderContext<OverridePixelOp>>(bitmap);//����Ԥ��Ⱦ�����ģ�����OverridePixelOp����������ʽ���ƣ�������ɫֱ�Ӹ���ԭ������ɫ 
        } else {//����͸��
            return makeAuto<PreviewRenderContext<OverPixelOp>>(bitmap);//����Ԥ��Ⱦ�����ģ�����OverPixelOp����������ʽ���ƣ���ԭ��ɫ�ϵ�������ɫ
        }
    }
}

void ParticleRenderer::render(const RenderParams& params, Statistics& stats, IRenderOutput& output) const {//������Ⱦ������Ⱦ�߼�
    MEASURE_SCOPE("ParticleRenderer::render");// ������Ⱦ�ܺ�ʱ

    Bitmap<Rgba> bitmap(params.camera->getSize());// ���������С������ȾĿ��λͼ
    AutoPtr<IRenderContext> context = getContext(params, bitmap);// ���ݲ�����ȡ��Ⱦ������

    // fill with the background color
    context->fill(params.background);// �ñ�����ɫ���λͼ

    if (grid > 0.f) {// ���ָ�����������
        drawGrid(*context, *params.camera, grid);//���Ƹ�������
    }

    struct {
        Vector r;//λ������
        Vector v;//��ɫ������
        bool used = false;//ʹ�ñ�־����¼�Ƿ�ѡ��
    } dir;//����dir�ṹ��������ʱ�洢�������ӵ�λ�ú��ٶ�,�Ӷ�ʵ�ֶ�ѡ�����ӵ����⴦��

    context->setColor(Rgba::black(), ColorFlag::LINE);// ��������Ĭ��������ɫΪ��ɫ

    shouldContinue = true;//����shouldContinue��־Ϊtrue����ʾ������Ⱦ
    // draw particles ��������

    const bool reverseOrder = dot(cached.cameraDir, params.camera->getFrame().row(2)) < 0._f;//ͨ�������жϱ������ӵ�˳���Ƿ���Ҫ��ת
    //ǰ��cameraDir������ӽǵķ��������������ǻ�ȡ�������ϵ�µ����Ϸ���������getFrame() ������������ϵ����row(2) ȡ����ĵ�3��,����ȡ���������ϵ�µ����Ϸ�������
    //������ߵĵ����������С��0,��ʾ�������������෴,Ҳ��������ӽǷ�ת��180�ȣ���ʱ�� reverseOrder ��Ϊ true��
    //�ڱ���������Ⱦʱ,��� reverseOrder Ϊ true,��ת���ӱ���˳�򣬴����һ�����ӿ�ʼ�������������Ա�֤�������ӻ�����Զ�����ӵ����棬�Ӷ�ʵ����ȷ�������ڵ���ϵ��

    for (Size k = 0; k < cached.positions.size(); ++k) {// ������������
        const Size i = reverseOrder ? cached.positions.size() - k - 1 : k; // �������ӱ���˳�򣨱������Ǹ�������ӽǣ�,������������
        if (!params.particles.renderGhosts && cached.idxs[i] == GHOST_INDEX) {// ��������������Ҳ���Ⱦ��������,������
            continue;
        }
        if (params.particles.selected && cached.idxs[i] == params.particles.selected.value()) {// �����ѡ�е�����
            // highlight the selected particle ����ѡ������
            context->setColor(Rgba::red(), ColorFlag::FILL);//���������ɫΪ��ɫ
            context->setColor(Rgba::white(), ColorFlag::LINE);//����������ɫΪ��ɫ

            if (!cached.vectors.empty()) {//��������е���ɫ���������ݲ�Ϊ��
                dir.used = true;//����ʹ�ñ�־Ϊtrue������ѡ��
                dir.v = cached.vectors[i];// �洢��ǰ���ӵ���ɫ������
                dir.r = cached.positions[i];// �洢��ǰ���ӵ�λ������
            }
        } else {//�������ѡ������
            Rgba color = cached.colors[i];// ��ȡ���ӵ���ɫ
            if (params.particles.grayScale) {// �����Ҫ�ҶȻ�
                color = Rgba(color.intensity());// ת��Ϊ�Ҷ���ɫ
            }
            context->setColor(color, ColorFlag::FILL | ColorFlag::LINE);// �������ӵ������ɫ��������ɫ
            if (cached.idxs[i] == GHOST_INDEX || cached.idxs[i] == ATTRACTOR_INDEX) {//������������ӻ�������
                context->setColor(Rgba::gray(0.7f), ColorFlag::LINE);// ����ǳ��ɫ����
            }
        }

        const Optional<ProjectedPoint> p = params.camera->project(cached.positions[i]);// ͶӰ���ӵ���Ļ�ռ�
        if (p) {//���ͶӰ�ɹ�
            float radius;//���Ӱ뾶
            if (cached.idxs[i] != ATTRACTOR_INDEX) {// �������������
                radius = p->radius * params.particles.scale;// ������������Ļ�ռ�İ뾶������������
            } else {//��������
                radius = p->radius;//������
            }
            const float size = min<float>(radius, context->size().x);// �����ڻ������ߴ�����
            context->drawCircle(p->coords, size);// ��������
        }
    }
    // after all particles are drawn, draw the velocity vector over �������������Ӻ󣬻����ٶ�ʸ��
    if (dir.used) {//����Ǳ�ѡ�е�����
        drawVector(*context, *params.camera, dir.r, dir.v, params.vectors.length);// ����ѡ�����ӵ��ٶ�����
    }

    renderOverlay(*context, params, stats);// ��Ⱦ���ǲ�

    // lastly black frame to draw on top of other stuff �����ƺ�ɫ�߿�,����ͼƬ�������ݵ��ϲ�
    const Pixel upper = bitmap.size() - Pixel(1, 1);// ����ͼ���С,��ȥ1��������Ϊ�߿���
    context->setColor(Rgba::black(), ColorFlag::LINE);// ���ñ߿���ɫΪ��ɫ
    //���ƾ��α߿��������
    context->drawLine(Coords(0, 0), Coords(upper.x, 0));
    context->drawLine(Coords(upper.x, 0), Coords(upper));
    context->drawLine(Coords(upper), Coords(0, upper.y));
    context->drawLine(Coords(0, upper.y), Coords(0, 0));

    output.update(bitmap, context->getLabels(), true);// �������
}

void ParticleRenderer::cancelRender() {//ȡ��������Ⱦ
    shouldContinue = false;//��shouldContinue����Ϊfalse����ʾȡ����Ⱦ
}

NAMESPACE_SPH_END
