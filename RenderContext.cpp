#include "objects/RenderContext.h"

NAMESPACE_SPH_BEGIN
//RenderContext渲染上下文：是渲染的前期准备接口类，内含一些公用的基础函数，来提供给大部分渲染器使用。
template <typename PixelOp>
void PreviewRenderContext<PixelOp>::setColor(const Rgba& color, const Flags<ColorFlag> flags) {//设置颜色
    //根据flags中的标志位,将对应的颜色成员变量设为传入的颜色。该函数允许同时设置多种颜色,通过flags进行控制。
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

//设置绘图属性,包括线条粗细、字体大小和位图填充颜色。
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
    bitmap.fill(color);//用于填充画布,以参数color作为背景色
}

template <typename PixelOp>
void PreviewRenderContext<PixelOp>::drawLine(Coords p1, Coords p2) {//绘制线条
    // 判断线条的斜率,如果横坐标变化较大则按x轴方向横向绘制,否则按y轴方向纵向绘制。
    if (abs(p2.x - p1.x) > abs(p2.y - p1.y)) {
        if (p1.x > p2.x) {
            std::swap(p1, p2);//交换起点和终点,保证绘制方向从小到大
        }
        const int x1 = int(floor(p1.x));//起点坐标向下取整
        const int x2 = int(ceil(p2.x));//终点坐标向上取整，从而获得整数型的起点和终点坐标
        const float y1 = p1.y;
        const float y2 = p2.y;
        for (int x = x1; x <= x2; ++x) {
            int y = int(y1 + (x - x1) * (y2 - y1) / (x2 - x1));//纵坐标为起点纵坐标加上横坐标变化量乘斜率
            drawSafe(Pixel(x, y), colors.line);//调用drawSafe在位图上绘制像素点，drawsafe不会绘制超出位图范围的点
        }
    } else {//同理，沿y轴方向绘制
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
void PreviewRenderContext<PixelOp>::drawCircle(const Coords center, const float radius) {//画圆
    if (center.x < -radius || center.x > bitmap.size().x + radius || center.y < -radius ||
        center.y > bitmap.size().y + radius) {//判断圆心x坐标是否小于负半径或者大于位图宽度+半径，y坐标是否小于负半径或者大于位图高度 + 半径
        //即判断该圆是否完全位于位图之外
        return;
    }
    const Pixel p(center);
    const int intRadius = min(int(radius), bitmap.size().x);//半径最大值为位图宽度
    if (p.x >= intRadius && p.x < bitmap.size().x - intRadius - 1 && p.y >= intRadius &&
        p.y < bitmap.size().y - intRadius - 1) {//判断圆心是否在足够绘制出完整圆形的区域内
        // can draw without checking 直接调用draw绘制，省去一次边界判断
        for (int y = -intRadius; y <= intRadius; ++y) {
            for (int x = -intRadius; x <= intRadius; ++x) {//遍历圆形的每个像素点
                const int rSqr = sqr(x) + sqr(y);//计算每点到圆心的距离
                if (rSqr <= sqr(radius - 1)) {//如果在内径内，即在圆内部，用填充颜色绘制
                    draw(p + Pixel(x, y), colors.fill);
                } else if (rSqr <= sqr(radius)) {//如果在外径内，即在圆的边界上，用线条颜色绘制
                    draw(p + Pixel(x, y), colors.line);
                }
            }
        }
    } else {//否则调用drawsafe绘制，多执行一次边界判断
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
void PreviewRenderContext<PixelOp>::drawTriangle(const Coords p1, const Coords p2, const Coords p3) {//绘制三角形：使用了扫描线填充算法
    // http://www-users.mat.uni.torun.pl/~wrona/3d_tutor/tri_fillers.html

    StaticArray<Coords, 3> p{ p1, p2, p3 };
    std::sort(p.begin(), p.end(), [](Coords p1, Coords p2) { return p1.y < p2.y; });//对p数组按y坐标从小到大排序
    SPH_ASSERT(p[0].y <= p[1].y && p[1].y <= p[2].y); // sanity check 保证p[0]的纵坐标最小,p[2]的纵坐标最大
    Coords a = p[0];//a为最低点，b为中间点，c为最高点
    Coords b = p[1];
    Coords c = p[2];
    a.y--;
    c.y++;//将a的y坐标减1，将c的y坐标加1
    //确保扫描线能够完整覆盖三角形的边界，从而实现完整的填充

    auto getDx = [](const Coords p1, const Coords p2) {//计算x关于y的斜率
        if (p2.y - p1.y > 0) {
            return float(p2.x - p1.x) / (p2.y - p1.y);
        } else {//两点纵坐标相同
            return 0.f;
        }
    };
    const float dx1 = getDx(a, b);
    const float dx2 = getDx(a, c);//dx2最大
    const float dx3 = getDx(b, c);

    auto doLine = [this](float x1, float x2, float y) {//绘制一条水平扫描线: x1、x2为扫描线段的左右端点横坐标，y为扫描线的纵坐标
        if (x1 > x2) {
            std::swap(x1, x2);//保证x1在左边,x2在右边
        }
        for (int x = int(floor(x1)); x <= int(ceil(x2)); ++x) {//x1向下取整，x2向上取整，遍历绘制
            drawSafe(Pixel(x, int(y)), colors.fill);
        }
    };

    Coords s = a, e = a;//s和e代表当前扫描线与两条边的交点
    //第一层循环扫描a到b的线段
    for (; s.y <= b.y; s.y++, e.y++, s.x += dx2, e.x += dx1) {//根据斜率dx计算对应的x坐标增量：相当于dx*1
        //这里默认了左端点在ac边上，右端点在ab边上。若不在，遍历时s+=dx2会使s>e,在doLine函数中会执行swap从而实现s和e的互换，确保他们位于正确的边上。
        //每次迭代s和e的y坐标增加1,即向上扫描一行
        doLine(s.x, e.x, s.y);//绘制扫描线
    }
    e = b;//扫描到b点后,将e置为b
    //第二层循环扫描b到c的线段
    for (; s.y <= c.y; s.y++, e.y++, s.x += dx2, e.x += dx3) {//默认左端点在ac边上，右端点在bc边上
        doLine(s.x, e.x, s.y);
    }
}

template <typename PixelOp>
void PreviewRenderContext<PixelOp>::drawBitmap(const Coords p, const Bitmap<Rgba>& subBitmap) {//绘制一个子位图，p表示子位图的绘制起点坐标，subBitmap表示需绘制的子绘图
    for (int y = 0; y < subBitmap.size().y; ++y) {
        for (int x = 0; x < subBitmap.size().x; ++x) {
            drawSafe(Pixel(x, y) + Pixel(p), subBitmap[Pixel(x, y)]);
        }
    }
}

template <typename PixelOp>
void PreviewRenderContext<PixelOp>::drawText(const Coords p, const Flags<TextAlign> align, const String& s) {//绘制文本，align参数指定文本的对齐方式
    labels.push(IRenderOutput::Label{ s, colors.text, fontSize, align, Pixel(p) });//不立即绘制,而是先将文本信息封装成一个Label对象。
    //接着将Label对象加入labels队列中缓存起来，实际的绘制会在后续的渲染阶段进行。
    //通过先缓存收集所有需要绘制的文本信息,然后在渲染时一次性绘制, 从而避免多次字符绘制带来的性能问题，以提升渲染性能。
}

//PreviewRenderContext模板类的模板实例化，获得了两种不同像素行为的渲染上下文
template class PreviewRenderContext<OverridePixelOp>;//覆盖式绘制：用新颜色直接覆盖原像素颜色 
template class PreviewRenderContext<OverPixelOp>;//叠加式绘制：在原颜色上叠加新颜色

void AntiAliasedRenderContext::drawCircle(const Coords center, const float radius) {//抗锯齿画圆
    if (center.x < -radius || center.x > bitmap.size().x + radius || center.y < -radius ||
        center.y > bitmap.size().y + radius) {//判断该圆是否完全位于位图之外
        return;
    }
    const Pixel p(center);
    if (radius <= 1.f) {//半径<=1时直接绘制圆心像素
        Rgba color = colors.fill;
        color.a() = sqr(radius);//alpha值即透明度，alpha值为0时完全透明，为1时完全不透明。
        //设置透明度为半径的平方，半径越短该点越透明。
        drawSafe(p, color);
    } else {
        const int r = int(std::ceil(radius)) + 1;// r取圆半径向上取整加1,保证覆盖圆形所有像素。
        for (int y = p.y - r; y <= p.y + r; ++y) {
            for (int x = p.x - r; x <= p.x + r; ++x) {//循环遍历以圆心为中心,半径向外扩展r的正方形区域。
                const float distSqr = sqr(x - center.x) + sqr(y - center.y);
                Rgba color = colors.fill;
                color.a() = clamp(radius - sqrt(distSqr), 0.f, 1.f);//设置alpha值为半径减去距离，同时使用clamp函数将其限制在0到1之间。
                //距离圆心越远的点，透明度越低—>越靠近边缘越透明，从而实现抗锯齿化的效果
                drawSafe(Pixel(x, y), color);
            }
        }
    }
}

void SmoothedRenderContext::drawCircle(const Coords center, const float radius) {//也是抗锯齿画圆，区别在于借助LutKernel核函数实现平滑化处理。
    //采用高斯核函数作为距离圆心的权重函数,可以获得更高质量的抗锯齿效果。
    if (center.x < -radius || center.x > bitmap.size().x + radius || center.y < -radius ||
        center.y > bitmap.size().y + radius) {
        return;
    }
    const Pixel p(center);
    const float maxRadius = radius * float(kernel.radius());//渲染时需要遍历超过圆形边缘一定范围的像素，以绘制高质量的抗锯齿圆形边缘。
    const float normalization = 1.f / float(kernel.valueImpl(0)); // sqr(25.f / particleScale); 计算了核函数的归一化系数。
    const int r = int(std::ceil(maxRadius)) + 1;//取半径为最大半径向上取整再加1，确保半径完全覆盖圆
    for (int y = p.y - r; y <= p.y + r; ++y) {
        for (int x = p.x - r; x <= p.x + r; ++x) {
            const float distSqr = sqr(x - center.x) + sqr(y - center.y);
            if (distSqr <= sqr(maxRadius + 1)) {//圆内
                Rgba color = colors.fill;

                const float alpha = float(kernel.valueImpl(distSqr / sqr(radius))) * normalization;//根据距离圆心的举例,计算核函数值。
                //核函数值经归一化到0-1范围，归一化后作为alpha通道值。从而实现距离圆心更远的像素,透明度更小,产生平滑效果。
                color.a() = clamp(alpha, 0.f, 1.f);//使用clamp函数再次确保alpha值在0-1范围内。
                drawSafe(Pixel(x, y), color);
            }
        }
    }
}

NAMESPACE_SPH_END
