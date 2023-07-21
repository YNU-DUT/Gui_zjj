#include "gui/objects/Palette.h"
#include "gui/objects/RenderContext.h"
#include "io/Path.h"
#include "objects/utility/Streams.h"
#include "post/Plot.h"
#include "post/Point.h"

NAMESPACE_SPH_BEGIN

float Palette::paletteToLinear(const float value) const {
    float palette;
    switch (scale) {// ���ݵ�ɫ��̶ȳ߼����ɫ��ֵ
    case PaletteScale::LINEAR:// ����ģʽ
        palette = value;//ֱ�ӷ�������ֵ
        break;
    case PaletteScale::LOGARITHMIC:// ����ģʽ
        // we allow calling this function with zero or negative value, it should simply map to the lowest
        // value on the palette ����������С�ڵ���0��ֵ���ô˺�������Ӧ�ü򵥵�ӳ�䵽��ɫ���ϵ����ֵ
        if (value < EPS) {//С��EPS������ֵӳ��Ϊ��Сֵ
            palette = -LARGE;
        } else {
            palette = float(log10(value));//��10Ϊ��ȡ����
        }
        break;
    case PaletteScale::HYBRID:// ���ģʽ
        if (value > 1.f) {//����ֵ����1ʱ��10Ϊ��ȡ������1
            palette = 1.f + float(log10(value));
        } else if (value < -1.f) {//����ֵС��-1ʱ��10Ϊ��ȡ������ȥ1
            palette = -1.f - float(log10(-value));
        } else {//����ֵ��-1��1֮��ʱֱ�ӷ���
            palette = value;
        }
        break;
    default:
        NOT_IMPLEMENTED;
    }
    SPH_ASSERT(isReal(palette), value);// ����palette��һ��ʵ��
    return palette;
}

float Palette::linearToPalette(const float value) const {
    switch (scale) {// ���ݵ�ɫ��̶ȳߴ����Կռ�ӳ�䵽��ɫ��ռ�
    case PaletteScale::LINEAR:// ����ģʽֱ�ӷ�������ֵ
        return value;
    case PaletteScale::LOGARITHMIC:// ����ģʽ:ȡ��10Ϊ������ֵ��ָ��
        return float(exp10(value));
    case PaletteScale::HYBRID:// ���ģʽ
        if (value > 1.f) {//����ֵ����1ʱ��10Ϊ��ȡָ����1
            return float(exp10(value - 1.f));
        } else if (value < -1.f) {//����ֵС��-1ʱ��10Ϊ��ȡָ����ȥ1
            return float(-exp10(-value - 1.f));
        } else {//����ֵ��-1��1֮��ʱֱ�ӷ���
            return value;
        }
    default:
        NOT_IMPLEMENTED; // in case new scale is added ��������¿̶�
    }
}

Palette::Palette(const Palette& other)// �������캯��
    : points(other.points.clone())// ��¡������ɫ��Ŀ��Ƶ㼯
    , range(other.range)// ���Ʒ�Χ
    , scale(other.scale) {}// ��������ģʽ

Palette& Palette::operator=(const Palette& other) {// ��ֵ���������
    points = copyable(other.points);//����������ɫ��Ŀ��Ƶ㼯
    range = other.range;// ���Ʒ�Χ
    scale = other.scale;// ��������ģʽ
    return *this;
}

Palette::Palette(Array<Point>&& controlPoints, const Interval& range, const PaletteScale scale)// ���캯��
    : points(std::move(controlPoints))// �ƶ����Ƶ㼯
    , range(range)// ���÷�Χ
    , scale(scale) {// ��������ģʽ
#ifdef SPH_DEBUG// ����ģʽ
    SPH_ASSERT(points.size() >= 2);// ���Կ��Ƶ㼯��������2����
    if (scale == PaletteScale::LOGARITHMIC) {//����ģʽ��
        SPH_ASSERT(range.lower() > 0.f);// ���Է�Χ���޴���0
    }
    SPH_ASSERT(std::all_of(points.begin(), points.end(), [](const Point& p) { // �������п��Ƶ��ֵ��0��1֮��
        return p.value >= 0 && p.value <= 1;
    }));
    // sanity check, points must be sorted �����Լ�飬���������
    SPH_ASSERT(std::is_sorted(points.begin(), points.end()));//���Կ��Ƶ��Ѿ�����
#endif
}

void Palette::addFixedPoint(const float value, const Rgba color) {//����¿��Ƶ�
    /// \todo store separately, do not move in setInterval! �ֿ��洢�������趨������ƶ�!
    points.push(Point{ this->rangeToRelative(value), color });//��㼯�в����¿��Ƶ�
    std::sort(points.begin(), points.end());//��points��������
}

const Array<Palette::Point>& Palette::getPoints() const {//��ȡ���Ƶ���Ϣ
    return points;
}

Interval Palette::getInterval() const {// ��ȡ��Χ
    SPH_ASSERT(points.size() >= 2);//����������2�����Ƶ�
    return range;
}

void Palette::setInterval(const Interval& newRange) {//���÷�Χ
    range = newRange;
}

PaletteScale Palette::getScale() const {//��ȡ����ģʽ
    return scale;
}

void Palette::setScale(const PaletteScale& newScale) {//��������ģʽ
    scale = newScale;
}

Rgba Palette::operator()(const float value) const {// �����()����:����ֵ���ض�Ӧ�ĵ�ɫ����ɫ
    const float x = rangeToRelative(value);// ��value�ľ���ֵӳ�䵽0-1���䣬��Ϊx
    SPH_ASSERT(points.size() >= 2);//����������2�����Ƶ�
    auto iter = std::lower_bound(points.begin(), points.end(), x, [](const Point& p, const float pos) { // ʹ��lower_bound���ҵ�һ����С��x�Ŀ��Ƶ����������Ϊiter
        return p.value < pos;
    });
    if (iter == points.begin()) {//iter���ڵ�һ�����Ƶ㣬Ҳ����xС�ڵ�һ�����Ƶ�
        return points.front().color;//���ص�һ�����Ƶ����ɫ
    } else if (iter == points.end()) {//iter�������һ�����Ƶ㣬Ҳ����x�������һ�����Ƶ�
        return points.back().color;//�������һ�����Ƶ����ɫ
    } else {// ��x�ڵ�һ�������һ�����Ƶ�֮��,�������Բ�ֵ�����Ӧ��ɫ
        const Rgba color2 = iter->color; // �õ� iter ָ��Ŀ��Ƶ���ɫ
        const Rgba color1 = (iter - 1)->color;// �õ� iter ǰһ�����Ƶ���ɫ
        const double pos2 = iter->value;// pos2ָ��iter��ָ��Ŀ��Ƶ��λ��
        const double pos1 = (iter - 1)->value;//pos1ָ��iter��ָ���ǰһ�����Ƶ��λ��
        const double f = (x - pos1) / (pos2 - pos1); // �����ֵ����
        return lerp(color1, color2, f);//ʹ��lerp�����������Ƶ���ɫ���в�ֵ,�Ӷ��õ�x��Ӧ����ɫ
        // lerp������������:
        // Rgba lerp(const Rgba& c1, const Rgba& c2, float f) {return c1 * (1 - f) + c2 * f;}
        // ������Ƶ���0.2��0.8λ��,��Ӧ��ɫΪred��blue:
        // ��x=0.5ʱ,f=(0.5-0.2)/(0.8-0.2)=0.5
        // �򷵻ص���ɫ��lerp(red, blue, 0.5)=purple

    }
}

Palette Palette::transform(Function<Rgba(const Rgba&)> func) const {// ��ɫ����ɫ�任
    Palette cloned = *this;// �����ǰ��ɫ��
    for (auto& point : cloned.points) {// �Ը�����ÿ�����Ƶ���ɫӦ�ô����func����
        point.color = func(point.color);
    }
    return cloned; // ���ظ���
}

Palette Palette::subsample(const Size pointCnt) const {// ��ɫ�����
    Array<Palette::Point> subsampled;// ���ѡȡ�Ŀ��Ƶ�
    for (Size i = 0; i < pointCnt; ++i) {
        const float x = float(i) / (pointCnt - 1);// ѡ���i�����Ƶ�,xΪ�ÿ��Ƶ���0-1�����λ��
        const float v = this->relativeToRange(x);// ��x��ӳ�䵽ֵ��Χ��
        subsampled.push(Palette::Point{ x, (*this)(v) });// ���ѡȡ�Ŀ��Ƶ�,��ɫͨ��ԭ��ɫ�����
    }
    return Palette(std::move(subsampled), this->getInterval(), this->getScale());// ����һ���µĵ�ɫ�壬Ϊԭ��ɫ��Ĳ��� 
}

float Palette::relativeToRange(const float value) const {// ��0-1�����value��ӳ�䵽��ɫ��ķ�Χ
    SPH_ASSERT(value >= 0.f && value <= 1.f);//����value��0-1֮��
    const float x1 = paletteToLinear(range.lower());// ��ɫ�巶Χ���������Կռ��е�λ��
    const float x2 = paletteToLinear(range.upper());// ��ɫ�巶Χ���������Կռ��е�λ��
    const float x = lerp(x1, x2, value);//ʹ��lerp����x1��x2���в�ֵ,�õ���Ӧvalue�����Կռ��е�λ��x
    return linearToPalette(x);// ��x��ӳ�䵽��ɫ��ķ�Χ 
}

float Palette::rangeToRelative(const float value) const {// ����ɫ���ֵ���ֵvalueӳ�䵽0 - 1����
    const float x1 = paletteToLinear(range.lower());// ��ɫ�巶Χ���������Կռ��е�λ��
    const float x2 = paletteToLinear(range.upper());// ��ɫ�巶Χ���������Կռ��е�λ��
    const float x = paletteToLinear(value);// value�����Կռ��е�λ��
    const float rel = (x - x1) / (x2 - x1);// relΪx��value��x1��x2֮������λ��
    SPH_ASSERT(isReal(rel), rel);//����rel��ʵ��
    return rel;
}

bool Palette::empty() const {//����ɫ���Ƿ�Ϊ��
    return points.empty();
}

Outcome Palette::loadFromStream(ITextInputStream& ifs) {// ���ı���ifs���ص�ɫ��
    try {
        Array<Rgba> colors;
        String line;
        while (ifs.readLine(line)) {
            line.replaceAll(",", " ");
            std::wstringstream ss(line.toUnicode());
            Rgba color;
            while (!ss.eof()) {
                ss >> color.r() >> color.g() >> color.b();
                color.a() = 1.f;
                colors.push(color);
            }
        }
        if (colors.size() < 2) {
            return makeFailed("No data loaded");
        }

        
        float from = points.front().value;
        float to = points.back().value;
        points.resize(colors.size());
        for (Size i = 0; i < points.size(); ++i) {
            points[i].color = colors[i];
            // yes, do not use linearToPalette, we want to map the palette in linear, not on quantities
            points[i].value = from + float(i) * (to - from) / (points.size() - 1);
        }
        return SUCCESS;
    } catch (const std::exception& e) {
        return makeFailed("Cannot load palette: {}", exceptionMessage(e));
    }
}

Outcome Palette::loadFromFile(const Path& path) {//���ļ�·�����ص�ɫ��
    FileTextInputStream ifs(path);
    return loadFromStream(ifs);
}

Outcome Palette::saveToStream(ITextOutputStream& ofs, const Size lineCnt) const {//�����ɫ�嵽�ı���ofs,ÿ��lineCnt����ɫ
    try {
        for (Size i = 0; i < lineCnt; ++i) {
            const float value = this->relativeToRange(float(i) / (lineCnt - 1));
            const Rgba color = this->operator()(value);
            ofs.write(format("{},{},{}", color.r(), color.g(), color.b()));
            if (i != lineCnt - 1) {
                ofs.write(L'\n');
            }
        }
        return SUCCESS;
    } catch (const std::exception& e) {
        return makeFailed("Cannot save palette: {}", exceptionMessage(e));
    }
}

Outcome Palette::saveToFile(const Path& path, const Size lineCnt) const {//�����ɫ�浽ָ���ļ�·��
    FileTextOutputStream ofs(path);
    return this->saveToStream(ofs, lineCnt);
}

void drawPalette(IRenderContext& context,//���Ƶ�ɫ��
    const Pixel origin,// ���
    const Pixel size,// ��С
    const Palette& palette,// ��ɫ��
    const Optional<Rgba>& lineColor) { // ��ѡ��������ɫ

    // draw palette
    for (int i = 0; i < size.x; ++i) {
        const float value = palette.relativeToRange(float(i) / (size.x - 1)); // ��ӳ���i�����ֵ
        context.setColor(palette(value), ColorFlag::LINE);// ���õ�i�������ɫ
        context.drawLine(Coords(origin.x + i, origin.y), Coords(origin.x + i, origin.y + size.y));//����һ����ֱ����
    }

    if (lineColor) {// ���ָ����������ɫ 
        // draw tics  ���ƿ̶� 
        const Interval interval = palette.getInterval();// ��ȡ��ɫ�巶Χ
        const PaletteScale scale = palette.getScale();// ��ȡ��ɫ��̶�

        Array<Float> tics;
        switch (scale) {// ���ݿ̶ȳ߼���̶�
        case PaletteScale::LINEAR://����ģʽ
            tics = getLinearTics(interval, 4);// ����4�����Կ̶�
            break;
        case PaletteScale::LOGARITHMIC: {//����ģʽ
            const Float lower = max(interval.lower(), 1.e-6_f);// ȡ���޺�1e-6֮��Ľϴ���
            const Float upper = interval.upper();// ����
            tics = getLogTics(Interval(lower, upper), 3, 5);// ��lower-upper�������3-5�������̶�
            break;
        }
        case PaletteScale::HYBRID: {//���ģʽ
            const Float lower = min(interval.lower(), -2._f);// ȡ���޺�-2֮��Ľ�С��
            const Float upper = max(interval.upper(), 2._f);// ȡ���޺�2֮��Ľϴ���
            tics = getHybridTics(Interval(lower, upper), 3);// ����չ����������3����Ͽ̶�
            break;
        }
        default:
            NOT_IMPLEMENTED;
        }
        context.setColor(lineColor.value(), ColorFlag::LINE | ColorFlag::TEXT);// ����������ɫ���ı���ɫ
        context.setFontSize(7); // ���������СΪ7
        for (Float tic : tics) {// �����̶ȼ���tics
            const float value = palette.rangeToRelative(float(tic));// ���̶�ֵticӳ�䵽0-1���� 
            const int i = int(value * size.x);// �����Ӧ��value��x����
            context.drawLine(Coords(origin.x + i, origin.y), Coords(origin.x + i, origin.y + 6));// �����Ͽ̶��� 
            context.drawLine(
                Coords(origin.x + i, origin.y + size.y - 6), Coords(origin.x + i, origin.y + size.y)); // �����¿̶��� 

            String text = toPrintableString(tic, 1, 10000);// ���̶�ֵ��ʽ��Ϊ�ַ���
            context.drawText(Coords(origin.x + i, origin.y + size.y + 15),// ���ƿ̶�ֵ�ı� 
                TextAlign::HORIZONTAL_CENTER | TextAlign::VERTICAL_CENTER,// ���в���
                text);
        }
    }
}


NAMESPACE_SPH_END
