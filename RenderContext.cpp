#include "objects/RenderContext.h"

NAMESPACE_SPH_BEGIN
//RenderContext��Ⱦ�����ģ�����Ⱦ��ǰ��׼���ӿ��࣬�ں�һЩ���õĻ������������ṩ���󲿷���Ⱦ��ʹ�á�
template <typename PixelOp>
void PreviewRenderContext<PixelOp>::setColor(const Rgba& color, const Flags<ColorFlag> flags) {//������ɫ
    //����flags�еı�־λ,����Ӧ����ɫ��Ա������Ϊ�������ɫ���ú�������ͬʱ���ö�����ɫ,ͨ��flags���п��ơ�
    if (flags.has(ColorFlag::LINE)) {
        colors.line = color;
    }
    if (flags.has(ColorFlag::FILL)) {
        colors.fill = color;
    }
    if (flags.has(ColorFlag::TEXT)) {
        colors.text = color;
    }
    // SPH_ASSERT(colors.line.alpha() == 1.f);
}

//���û�ͼ����,����������ϸ�������С��λͼ�����ɫ��
template <typename PixelOp>
void PreviewRenderContext<PixelOp>::setThickness(const float newThickness) {
    thickness = newThickness;
}

template <typename PixelOp>
void PreviewRenderContext<PixelOp>::setFontSize(const int newFontSize) {
    fontSize = newFontSize;
}

template <typename PixelOp>
void PreviewRenderContext<PixelOp>::fill(const Rgba& color) {
    bitmap.fill(color);//������仭��,�Բ���color��Ϊ����ɫ
}

template <typename PixelOp>
void PreviewRenderContext<PixelOp>::drawLine(Coords p1, Coords p2) {//��������
    // �ж�������б��,���������仯�ϴ���x�᷽��������,����y�᷽��������ơ�
    if (abs(p2.x - p1.x) > abs(p2.y - p1.y)) {
        if (p1.x > p2.x) {
            std::swap(p1, p2);//���������յ�,��֤���Ʒ����С����
        }
        const int x1 = int(floor(p1.x));//�����������ȡ��
        const int x2 = int(ceil(p2.x));//�յ���������ȡ�����Ӷ���������͵������յ�����
        const float y1 = p1.y;
        const float y2 = p2.y;
        for (int x = x1; x <= x2; ++x) {
            int y = int(y1 + (x - x1) * (y2 - y1) / (x2 - x1));//������Ϊ�����������Ϻ�����仯����б��
            drawSafe(Pixel(x, y), colors.line);//����drawSafe��λͼ�ϻ������ص㣬drawsafe������Ƴ���λͼ��Χ�ĵ�
        }
    } else {//ͬ����y�᷽�����
        if (p1.y > p2.y) {
            std::swap(p1, p2);
        }
        const int y1 = int(floor(p1.y));
        const int y2 = int(ceil(p2.y));
        const float x1 = p1.x;
        const float x2 = p2.x;
        for (int y = y1; y <= y2; ++y) {
            int x = int(x1 + (y - y1) * (x2 - x1) / (y2 - y1));
            drawSafe(Pixel(x, y), colors.line);
        }
    }
}

template <typename PixelOp>
void PreviewRenderContext<PixelOp>::drawCircle(const Coords center, const float radius) {//��Բ
    if (center.x < -radius || center.x > bitmap.size().x + radius || center.y < -radius ||
        center.y > bitmap.size().y + radius) {//�ж�Բ��x�����Ƿ�С�ڸ��뾶���ߴ���λͼ���+�뾶��y�����Ƿ�С�ڸ��뾶���ߴ���λͼ�߶� + �뾶
        //���жϸ�Բ�Ƿ���ȫλ��λͼ֮��
        return;
    }
    const Pixel p(center);
    const int intRadius = min(int(radius), bitmap.size().x);//�뾶���ֵΪλͼ���
    if (p.x >= intRadius && p.x < bitmap.size().x - intRadius - 1 && p.y >= intRadius &&
        p.y < bitmap.size().y - intRadius - 1) {//�ж�Բ���Ƿ����㹻���Ƴ�����Բ�ε�������
        // can draw without checking ֱ�ӵ���draw���ƣ�ʡȥһ�α߽��ж�
        for (int y = -intRadius; y <= intRadius; ++y) {
            for (int x = -intRadius; x <= intRadius; ++x) {//����Բ�ε�ÿ�����ص�
                const int rSqr = sqr(x) + sqr(y);//����ÿ�㵽Բ�ĵľ���
                if (rSqr <= sqr(radius - 1)) {//������ھ��ڣ�����Բ�ڲ����������ɫ����
                    draw(p + Pixel(x, y), colors.fill);
                } else if (rSqr <= sqr(radius)) {//������⾶�ڣ�����Բ�ı߽��ϣ���������ɫ����
                    draw(p + Pixel(x, y), colors.line);
                }
            }
        }
    } else {//�������drawsafe���ƣ���ִ��һ�α߽��ж�
        for (int y = -intRadius; y <= intRadius; ++y) {
            for (int x = -intRadius; x <= intRadius; ++x) {
                const int rSqr = sqr(x) + sqr(y);
                if (rSqr <= sqr(radius - 1)) {
                    drawSafe(p + Pixel(x, y), colors.fill);
                } else if (rSqr <= sqr(radius)) {
                    drawSafe(p + Pixel(x, y), colors.line);
                }
            }
        }
    }
}

template <typename PixelOp>
void PreviewRenderContext<PixelOp>::drawTriangle(const Coords p1, const Coords p2, const Coords p3) {//���������Σ�ʹ����ɨ��������㷨
    // http://www-users.mat.uni.torun.pl/~wrona/3d_tutor/tri_fillers.html

    StaticArray<Coords, 3> p{ p1, p2, p3 };
    std::sort(p.begin(), p.end(), [](Coords p1, Coords p2) { return p1.y < p2.y; });//��p���鰴y�����С��������
    SPH_ASSERT(p[0].y <= p[1].y && p[1].y <= p[2].y); // sanity check ��֤p[0]����������С,p[2]�����������
    Coords a = p[0];//aΪ��͵㣬bΪ�м�㣬cΪ��ߵ�
    Coords b = p[1];
    Coords c = p[2];
    a.y--;
    c.y++;//��a��y�����1����c��y�����1
    //ȷ��ɨ�����ܹ��������������εı߽磬�Ӷ�ʵ�����������

    auto getDx = [](const Coords p1, const Coords p2) {//����x����y��б��
        if (p2.y - p1.y > 0) {
            return float(p2.x - p1.x) / (p2.y - p1.y);
        } else {//������������ͬ
            return 0.f;
        }
    };
    const float dx1 = getDx(a, b);
    const float dx2 = getDx(a, c);//dx2���
    const float dx3 = getDx(b, c);

    auto doLine = [this](float x1, float x2, float y) {//����һ��ˮƽɨ����: x1��x2Ϊɨ���߶ε����Ҷ˵�����꣬yΪɨ���ߵ�������
        if (x1 > x2) {
            std::swap(x1, x2);//��֤x1�����,x2���ұ�
        }
        for (int x = int(floor(x1)); x <= int(ceil(x2)); ++x) {//x1����ȡ����x2����ȡ������������
            drawSafe(Pixel(x, int(y)), colors.fill);
        }
    };

    Coords s = a, e = a;//s��e����ǰɨ�����������ߵĽ���
    //��һ��ѭ��ɨ��a��b���߶�
    for (; s.y <= b.y; s.y++, e.y++, s.x += dx2, e.x += dx1) {//����б��dx�����Ӧ��x�����������൱��dx*1
        //����Ĭ������˵���ac���ϣ��Ҷ˵���ab���ϡ������ڣ�����ʱs+=dx2��ʹs>e,��doLine�����л�ִ��swap�Ӷ�ʵ��s��e�Ļ�����ȷ������λ����ȷ�ı��ϡ�
        //ÿ�ε���s��e��y��������1,������ɨ��һ��
        doLine(s.x, e.x, s.y);//����ɨ����
    }
    e = b;//ɨ�赽b���,��e��Ϊb
    //�ڶ���ѭ��ɨ��b��c���߶�
    for (; s.y <= c.y; s.y++, e.y++, s.x += dx2, e.x += dx3) {//Ĭ����˵���ac���ϣ��Ҷ˵���bc����
        doLine(s.x, e.x, s.y);
    }
}

template <typename PixelOp>
void PreviewRenderContext<PixelOp>::drawBitmap(const Coords p, const Bitmap<Rgba>& subBitmap) {//����һ����λͼ��p��ʾ��λͼ�Ļ���������꣬subBitmap��ʾ����Ƶ��ӻ�ͼ
    for (int y = 0; y < subBitmap.size().y; ++y) {
        for (int x = 0; x < subBitmap.size().x; ++x) {
            drawSafe(Pixel(x, y) + Pixel(p), subBitmap[Pixel(x, y)]);
        }
    }
}

template <typename PixelOp>
void PreviewRenderContext<PixelOp>::drawText(const Coords p, const Flags<TextAlign> align, const String& s) {//�����ı���align����ָ���ı��Ķ��뷽ʽ
    labels.push(IRenderOutput::Label{ s, colors.text, fontSize, align, Pixel(p) });//����������,�����Ƚ��ı���Ϣ��װ��һ��Label����
    //���Ž�Label�������labels�����л���������ʵ�ʵĻ��ƻ��ں�������Ⱦ�׶ν��С�
    //ͨ���Ȼ����ռ�������Ҫ���Ƶ��ı���Ϣ,Ȼ������Ⱦʱһ���Ի���, �Ӷ��������ַ����ƴ������������⣬��������Ⱦ���ܡ�
}

//PreviewRenderContextģ�����ģ��ʵ��������������ֲ�ͬ������Ϊ����Ⱦ������
template class PreviewRenderContext<OverridePixelOp>;//����ʽ���ƣ�������ɫֱ�Ӹ���ԭ������ɫ 
template class PreviewRenderContext<OverPixelOp>;//����ʽ���ƣ���ԭ��ɫ�ϵ�������ɫ

void AntiAliasedRenderContext::drawCircle(const Coords center, const float radius) {//����ݻ�Բ
    if (center.x < -radius || center.x > bitmap.size().x + radius || center.y < -radius ||
        center.y > bitmap.size().y + radius) {//�жϸ�Բ�Ƿ���ȫλ��λͼ֮��
        return;
    }
    const Pixel p(center);
    if (radius <= 1.f) {//�뾶<=1ʱֱ�ӻ���Բ������
        Rgba color = colors.fill;
        color.a() = sqr(radius);//alphaֵ��͸���ȣ�alphaֵΪ0ʱ��ȫ͸����Ϊ1ʱ��ȫ��͸����
        //����͸����Ϊ�뾶��ƽ�����뾶Խ�̸õ�Խ͸����
        drawSafe(p, color);
    } else {
        const int r = int(std::ceil(radius)) + 1;// rȡԲ�뾶����ȡ����1,��֤����Բ���������ء�
        for (int y = p.y - r; y <= p.y + r; ++y) {
            for (int x = p.x - r; x <= p.x + r; ++x) {//ѭ��������Բ��Ϊ����,�뾶������չr������������
                const float distSqr = sqr(x - center.x) + sqr(y - center.y);
                Rgba color = colors.fill;
                color.a() = clamp(radius - sqrt(distSqr), 0.f, 1.f);//����alphaֵΪ�뾶��ȥ���룬ͬʱʹ��clamp��������������0��1֮�䡣
                //����Բ��ԽԶ�ĵ㣬͸����Խ�͡�>Խ������ԵԽ͸�����Ӷ�ʵ�ֿ���ݻ���Ч��
                drawSafe(Pixel(x, y), color);
            }
        }
    }
}

void SmoothedRenderContext::drawCircle(const Coords center, const float radius) {//Ҳ�ǿ���ݻ�Բ���������ڽ���LutKernel�˺���ʵ��ƽ��������
    //���ø�˹�˺�����Ϊ����Բ�ĵ�Ȩ�غ���,���Ի�ø��������Ŀ����Ч����
    if (center.x < -radius || center.x > bitmap.size().x + radius || center.y < -radius ||
        center.y > bitmap.size().y + radius) {
        return;
    }
    const Pixel p(center);
    const float maxRadius = radius * float(kernel.radius());//��Ⱦʱ��Ҫ��������Բ�α�Եһ����Χ�����أ��Ի��Ƹ������Ŀ����Բ�α�Ե��
    const float normalization = 1.f / float(kernel.valueImpl(0)); // sqr(25.f / particleScale); �����˺˺����Ĺ�һ��ϵ����
    const int r = int(std::ceil(maxRadius)) + 1;//ȡ�뾶Ϊ���뾶����ȡ���ټ�1��ȷ���뾶��ȫ����Բ
    for (int y = p.y - r; y <= p.y + r; ++y) {
        for (int x = p.x - r; x <= p.x + r; ++x) {
            const float distSqr = sqr(x - center.x) + sqr(y - center.y);
            if (distSqr <= sqr(maxRadius + 1)) {//Բ��
                Rgba color = colors.fill;

                const float alpha = float(kernel.valueImpl(distSqr / sqr(radius))) * normalization;//���ݾ���Բ�ĵľ���,����˺���ֵ��
                //�˺���ֵ����һ����0-1��Χ����һ������Ϊalphaͨ��ֵ���Ӷ�ʵ�־���Բ�ĸ�Զ������,͸���ȸ�С,����ƽ��Ч����
                color.a() = clamp(alpha, 0.f, 1.f);//ʹ��clamp�����ٴ�ȷ��alphaֵ��0-1��Χ�ڡ�
                drawSafe(Pixel(x, y), color);
            }
        }
    }
}

NAMESPACE_SPH_END
