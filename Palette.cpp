#include "gui/objects/Palette.h"
#include "gui/objects/RenderContext.h"
#include "io/Path.h"
#include "objects/utility/Streams.h"
#include "post/Plot.h"
#include "post/Point.h"

NAMESPACE_SPH_BEGIN

float Palette::paletteToLinear(const float value) const {
    float palette;
    switch (scale) {// 根据调色版刻度尺计算调色板值
    case PaletteScale::LINEAR:// 线性模式
        palette = value;//直接返回输入值
        break;
    case PaletteScale::LOGARITHMIC:// 对数模式
        // we allow calling this function with zero or negative value, it should simply map to the lowest
        // value on the palette 我们允许以小于等于0的值调用此函数，它应该简单地映射到调色板上的最低值
        if (value < EPS) {//小于EPS的输入值映射为最小值
            palette = -LARGE;
        } else {
            palette = float(log10(value));//以10为底取对数
        }
        break;
    case PaletteScale::HYBRID:// 混合模式
        if (value > 1.f) {//输入值大于1时以10为底取对数加1
            palette = 1.f + float(log10(value));
        } else if (value < -1.f) {//输入值小于-1时以10为底取对数减去1
            palette = -1.f - float(log10(-value));
        } else {//输入值在-1到1之间时直接返回
            palette = value;
        }
        break;
    default:
        NOT_IMPLEMENTED;
    }
    SPH_ASSERT(isReal(palette), value);// 断言palette是一个实数
    return palette;
}

float Palette::linearToPalette(const float value) const {
    switch (scale) {// 根据调色版刻度尺从线性空间映射到调色板空间
    case PaletteScale::LINEAR:// 线性模式直接返回输入值
        return value;
    case PaletteScale::LOGARITHMIC:// 对数模式:取以10为底输入值的指数
        return float(exp10(value));
    case PaletteScale::HYBRID:// 混合模式
        if (value > 1.f) {//输入值大于1时以10为底取指数加1
            return float(exp10(value - 1.f));
        } else if (value < -1.f) {//输入值小于-1时以10为底取指数减去1
            return float(-exp10(-value - 1.f));
        } else {//输入值在-1到1之间时直接返回
            return value;
        }
    default:
        NOT_IMPLEMENTED; // in case new scale is added 可以添加新刻度
    }
}

Palette::Palette(const Palette& other)// 拷贝构造函数
    : points(other.points.clone())// 克隆其他调色板的控制点集
    , range(other.range)// 复制范围
    , scale(other.scale) {}// 复制缩放模式

Palette& Palette::operator=(const Palette& other) {// 赋值运算符重载
    points = copyable(other.points);//拷贝其他调色板的控制点集
    range = other.range;// 复制范围
    scale = other.scale;// 复制缩放模式
    return *this;
}

Palette::Palette(Array<Point>&& controlPoints, const Interval& range, const PaletteScale scale)// 构造函数
    : points(std::move(controlPoints))// 移动控制点集
    , range(range)// 设置范围
    , scale(scale) {// 设置缩放模式
#ifdef SPH_DEBUG// 调试模式
    SPH_ASSERT(points.size() >= 2);// 断言控制点集合至少有2个点
    if (scale == PaletteScale::LOGARITHMIC) {//对数模式下
        SPH_ASSERT(range.lower() > 0.f);// 断言范围下限大于0
    }
    SPH_ASSERT(std::all_of(points.begin(), points.end(), [](const Point& p) { // 断言所有控制点的值在0到1之间
        return p.value >= 0 && p.value <= 1;
    }));
    // sanity check, points must be sorted 完整性检查，点必须排序
    SPH_ASSERT(std::is_sorted(points.begin(), points.end()));//断言控制点已经排序
#endif
}

void Palette::addFixedPoint(const float value, const Rgba color) {//添加新控制点
    /// \todo store separately, do not move in setInterval! 分开存储，不在设定间隔内移动!
    points.push(Point{ this->rangeToRelative(value), color });//向点集中插入新控制点
    std::sort(points.begin(), points.end());//对points进行排序
}

const Array<Palette::Point>& Palette::getPoints() const {//获取控制点信息
    return points;
}

Interval Palette::getInterval() const {// 获取范围
    SPH_ASSERT(points.size() >= 2);//断言至少有2个控制点
    return range;
}

void Palette::setInterval(const Interval& newRange) {//设置范围
    range = newRange;
}

PaletteScale Palette::getScale() const {//获取缩放模式
    return scale;
}

void Palette::setScale(const PaletteScale& newScale) {//设置缩放模式
    scale = newScale;
}

Rgba Palette::operator()(const float value) const {// 运算符()重载:根据值返回对应的调色板颜色
    const float x = rangeToRelative(value);// 将value的绝对值映射到0-1区间，记为x
    SPH_ASSERT(points.size() >= 2);//断言至少有2个控制点
    auto iter = std::lower_bound(points.begin(), points.end(), x, [](const Point& p, const float pos) { // 使用lower_bound查找第一个不小于x的控制点迭代器，记为iter
        return p.value < pos;
    });
    if (iter == points.begin()) {//iter等于第一个控制点，也就是x小于第一个控制点
        return points.front().color;//返回第一个控制点的颜色
    } else if (iter == points.end()) {//iter等于最后一个控制点，也就是x大于最后一个控制点
        return points.back().color;//返回最后一个控制点的颜色
    } else {// 若x在第一个和最后一个控制点之间,进行线性插值计算对应颜色
        const Rgba color2 = iter->color; // 得到 iter 指向的控制点颜色
        const Rgba color1 = (iter - 1)->color;// 得到 iter 前一个控制点颜色
        const double pos2 = iter->value;// pos2指向iter所指向的控制点的位置
        const double pos1 = (iter - 1)->value;//pos1指向iter所指向的前一个控制点的位置
        const double f = (x - pos1) / (pos2 - pos1); // 计算插值因子
        return lerp(color1, color2, f);//使用lerp函数对两控制点颜色进行插值,从而得到x对应的颜色
        // lerp函数定义如下:
        // Rgba lerp(const Rgba& c1, const Rgba& c2, float f) {return c1 * (1 - f) + c2 * f;}
        // 例如控制点在0.2和0.8位置,对应颜色为red和blue:
        // 当x=0.5时,f=(0.5-0.2)/(0.8-0.2)=0.5
        // 则返回的颜色是lerp(red, blue, 0.5)=purple

    }
}

Palette Palette::transform(Function<Rgba(const Rgba&)> func) const {// 调色板颜色变换
    Palette cloned = *this;// 深拷贝当前调色板
    for (auto& point : cloned.points) {// 对副本的每个控制点颜色应用传入的func函数
        point.color = func(point.color);
    }
    return cloned; // 返回副本
}

Palette Palette::subsample(const Size pointCnt) const {// 调色板采样
    Array<Palette::Point> subsampled;// 存放选取的控制点
    for (Size i = 0; i < pointCnt; ++i) {
        const float x = float(i) / (pointCnt - 1);// 选择第i个控制点,x为该控制点在0-1区间的位置
        const float v = this->relativeToRange(x);// 将x反映射到值域范围内
        subsampled.push(Palette::Point{ x, (*this)(v) });// 添加选取的控制点,颜色通过原调色板计算
    }
    return Palette(std::move(subsampled), this->getInterval(), this->getScale());// 返回一个新的调色板，为原调色板的采样 
}

float Palette::relativeToRange(const float value) const {// 将0-1区间的value反映射到调色板的范围
    SPH_ASSERT(value >= 0.f && value <= 1.f);//断言value在0-1之间
    const float x1 = paletteToLinear(range.lower());// 调色板范围下限在线性空间中的位置
    const float x2 = paletteToLinear(range.upper());// 调色板范围上限在线性空间中的位置
    const float x = lerp(x1, x2, value);//使用lerp函数x1、x2进行插值,得到对应value在线性空间中的位置x
    return linearToPalette(x);// 将x反映射到调色板的范围 
}

float Palette::rangeToRelative(const float value) const {// 将调色板的值域的值value映射到0 - 1区间
    const float x1 = paletteToLinear(range.lower());// 调色板范围下限在线性空间中的位置
    const float x2 = paletteToLinear(range.upper());// 调色板范围上限在线性空间中的位置
    const float x = paletteToLinear(value);// value在线性空间中的位置
    const float rel = (x - x1) / (x2 - x1);// rel为x即value在x1和x2之间的相对位置
    SPH_ASSERT(isReal(rel), rel);//断言rel是实数
    return rel;
}

bool Palette::empty() const {//检查调色板是否为空
    return points.empty();
}

Outcome Palette::loadFromStream(ITextInputStream& ifs) {// 从文本流ifs加载调色板
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

Outcome Palette::loadFromFile(const Path& path) {//从文件路径加载调色板
    FileTextInputStream ifs(path);
    return loadFromStream(ifs);
}

Outcome Palette::saveToStream(ITextOutputStream& ofs, const Size lineCnt) const {//保存调色板到文本流ofs,每行lineCnt个颜色
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

Outcome Palette::saveToFile(const Path& path, const Size lineCnt) const {//保存调色版到指定文件路径
    FileTextOutputStream ofs(path);
    return this->saveToStream(ofs, lineCnt);
}

void drawPalette(IRenderContext& context,//绘制调色板
    const Pixel origin,// 起点
    const Pixel size,// 大小
    const Palette& palette,// 调色板
    const Optional<Rgba>& lineColor) { // 可选的线条颜色

    // draw palette
    for (int i = 0; i < size.x; ++i) {
        const float value = palette.relativeToRange(float(i) / (size.x - 1)); // 反映射第i个点的值
        context.setColor(palette(value), ColorFlag::LINE);// 设置第i个点的颜色
        context.drawLine(Coords(origin.x + i, origin.y), Coords(origin.x + i, origin.y + size.y));//绘制一条垂直线条
    }

    if (lineColor) {// 如果指定了线条颜色 
        // draw tics  绘制刻度 
        const Interval interval = palette.getInterval();// 获取调色板范围
        const PaletteScale scale = palette.getScale();// 获取调色板刻度

        Array<Float> tics;
        switch (scale) {// 根据刻度尺计算刻度
        case PaletteScale::LINEAR://线性模式
            tics = getLinearTics(interval, 4);// 计算4个线性刻度
            break;
        case PaletteScale::LOGARITHMIC: {//对数模式
            const Float lower = max(interval.lower(), 1.e-6_f);// 取下限和1e-6之间的较大者
            const Float upper = interval.upper();// 上限
            tics = getLogTics(Interval(lower, upper), 3, 5);// 在lower-upper区间计算3-5个对数刻度
            break;
        }
        case PaletteScale::HYBRID: {//混合模式
            const Float lower = min(interval.lower(), -2._f);// 取下限和-2之间的较小者
            const Float upper = max(interval.upper(), 2._f);// 取上限和2之间的较大者
            tics = getHybridTics(Interval(lower, upper), 3);// 在扩展后的区间计算3个混合刻度
            break;
        }
        default:
            NOT_IMPLEMENTED;
        }
        context.setColor(lineColor.value(), ColorFlag::LINE | ColorFlag::TEXT);// 设置线条颜色和文本颜色
        context.setFontSize(7); // 设置字体大小为7
        for (Float tic : tics) {// 遍历刻度集合tics
            const float value = palette.rangeToRelative(float(tic));// 将刻度值tic映射到0-1区间 
            const int i = int(value * size.x);// 计算对应于value的x坐标
            context.drawLine(Coords(origin.x + i, origin.y), Coords(origin.x + i, origin.y + 6));// 绘制上刻度线 
            context.drawLine(
                Coords(origin.x + i, origin.y + size.y - 6), Coords(origin.x + i, origin.y + size.y)); // 绘制下刻度线 

            String text = toPrintableString(tic, 1, 10000);// 将刻度值格式化为字符串
            context.drawText(Coords(origin.x + i, origin.y + size.y + 15),// 绘制刻度值文本 
                TextAlign::HORIZONTAL_CENTER | TextAlign::VERTICAL_CENTER,// 居中布局
                text);
        }
    }
}


NAMESPACE_SPH_END
