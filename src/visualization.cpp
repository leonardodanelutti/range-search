#include "visualization.h"
#include <random>
#include <iostream>

RangeSearchWidget::RangeSearchWidget(QWidget* parent)
    : QWidget(parent)
    , dragging_(false)
    , resizing_(false)
    , showTreePartitions_(false)
    , pointCount_(20)
    , scale_(1.0)
    , offset_(0, 0)
{
    setMinimumSize(800, 600);
    setMouseTracking(true);
    
    // Initialize with a default query region
    queryRegion_ = QRectF(2.0, 2.0, 4.0, 3.0);
    
    // Generate some initial points
    generateRandomPoints();
}

void RangeSearchWidget::paintEvent(QPaintEvent* event) {
    // Init
    QPainter painter(this);
    painter.fillRect(rect(), Qt::white);
    painter.save();
    
    // Calculate scale to fit world coordinates into widget
    double widgetWidth = width() - 2 * MARGIN;
    double widgetHeight = height() - 2 * MARGIN;
    scale_ = std::min(widgetWidth / WORLD_SIZE, widgetHeight / WORLD_SIZE);
    
    offset_ = QPointF(MARGIN, MARGIN);
    
    // Draw coordinate axes and boundary
    painter.setPen(QPen(Qt::lightGray, 1));
    QPoint origin = worldToScreen(QPointF(0, 0));
    QPoint xEnd = worldToScreen(QPointF(WORLD_SIZE, 0));
    QPoint yEnd = worldToScreen(QPointF(0, WORLD_SIZE));
    QPoint topRight = worldToScreen(QPointF(WORLD_SIZE, WORLD_SIZE));
    QPoint topLeft = worldToScreen(QPointF(0, WORLD_SIZE));
    QPoint bottomRight = worldToScreen(QPointF(WORLD_SIZE, 0));
    
    painter.drawLine(origin, xEnd);
    painter.drawLine(origin, yEnd);
    painter.drawLine(topLeft, topRight);
    painter.drawLine(bottomRight, topRight);
    
    // Draw KD-tree partitions if enabled
    if (showTreePartitions_ && kdTree_ && !kdTree_->empty()) {
        drawKDTreePartitions(painter);
    }
    
    // Draw all points
    drawPoints(painter);
    
    // Draw query region
    drawQueryRegion(painter);
    
    // Draw result points (highlighted)
    drawResultPoints(painter);
    
    painter.restore();
}

void RangeSearchWidget::drawPoints(QPainter& painter) {
    painter.setPen(QPen(Qt::black, 2));
    painter.setBrush(Qt::black);
    
    for (const auto& point : points_) {
        double x = CGAL::to_double(point.cartesian(0));
        double y = CGAL::to_double(point.cartesian(1));
        QPoint screenPos = worldToScreen(QPointF(x, y));
        painter.drawEllipse(screenPos, 4, 4);
    }
}

void RangeSearchWidget::drawQueryRegion(QPainter& painter) {
    painter.setPen(QPen(Qt::red, 2));
    painter.setBrush(QBrush(Qt::red, Qt::Dense7Pattern));
    
    QPoint topLeft = worldToScreen(queryRegion_.topLeft());
    QPoint bottomRight = worldToScreen(queryRegion_.bottomRight());
    QRect screenRect(topLeft, bottomRight);
    
    painter.drawRect(screenRect);
    
    // Draw resize handle
    painter.setBrush(Qt::red);
    painter.drawEllipse(bottomRight.x() - 4, bottomRight.y() - 4, 8, 8);
}

void RangeSearchWidget::drawResultPoints(QPainter& painter) {
    painter.setPen(QPen(Qt::green, 3));
    painter.setBrush(Qt::green);
    
    for (const auto& point : searchResults_) {
        double x = CGAL::to_double(point.cartesian(0));
        double y = CGAL::to_double(point.cartesian(1));
        QPoint screenPos = worldToScreen(QPointF(x, y));
        painter.drawEllipse(screenPos, 6, 6);
    }
}

void RangeSearchWidget::drawKDTreePartitions(QPainter& painter) {
    if (!kdTree_ || kdTree_->empty()) return;
    
    painter.setPen(QPen(Qt::blue, 1, Qt::DashLine));
    drawTreePartitions(painter, kdTree_->getRoot(), 0, 0, WORLD_SIZE, WORLD_SIZE, 0);
}

void RangeSearchWidget::drawTreePartitions(QPainter& painter, KDTree<2>::Node* node,
                                           double minX, double minY, double maxX, double maxY, int depth) const {
    if (node->is_leaf()) return;
    
    KDTree<2>::InternalNode* internal = static_cast<KDTree<2>::InternalNode*>(node);
    int dimension = internal->split_dimension;
    double splitValue = internal->split_coordinate;
    
    if (dimension == 0) { // X
        if (splitValue >= minX && splitValue <= maxX) {
            QPoint start = worldToScreen(QPointF(splitValue, minY));
            QPoint end = worldToScreen(QPointF(splitValue, maxY));
            painter.drawLine(start, end);
            
            if (internal->left) {
                drawTreePartitions(painter, internal->left.get(), minX, minY, splitValue, maxY, depth + 1);
            }
            if (internal->right) {
                drawTreePartitions(painter, internal->right.get(), splitValue, minY, maxX, maxY, depth + 1);
            }
        }
    } else { // Y
        if (splitValue >= minY && splitValue <= maxY) {
            QPoint start = worldToScreen(QPointF(minX, splitValue));
            QPoint end = worldToScreen(QPointF(maxX, splitValue));
            painter.drawLine(start, end);
            
            if (internal->left) {
                drawTreePartitions(painter, internal->left.get(), minX, minY, maxX, splitValue, depth + 1);
            }
            if (internal->right) {
                drawTreePartitions(painter, internal->right.get(), minX, splitValue, maxX, maxY, depth + 1);
            }
        }
    }
}

QPointF RangeSearchWidget::screenToWorld(const QPoint& screen) const {
    double x = (screen.x() - offset_.x()) / scale_;
    double y = (screen.y() - offset_.y()) / scale_;
    return QPointF(x, y);
}

QPoint RangeSearchWidget::worldToScreen(const QPointF& world) const {
    int x = static_cast<int>(world.x() * scale_ + offset_.x());
    int y = static_cast<int>(world.y() * scale_ + offset_.y());
    return QPoint(x, y);
}

void RangeSearchWidget::mousePressEvent(QMouseEvent* event) {
    if (event->button() != Qt::LeftButton) return;
    
    QPointF worldPos = screenToWorld(event->pos());
    QPoint bottomRightScreen = worldToScreen(queryRegion_.bottomRight());
    
    // Check if clicking on resize handle
    if ((event->pos() - bottomRightScreen).manhattanLength() < 10) {
        resizing_ = true;
        dragStart_ = worldPos;
        originalTopLeft_ = queryRegion_.topLeft();
        originalSize_ = queryRegion_.size();
    }
    // Check if clicking inside query region for moving
    else if (queryRegion_.contains(worldPos)) {
        dragging_ = true;
        dragStart_ = worldPos;
        originalTopLeft_ = queryRegion_.topLeft();
    }
}

void RangeSearchWidget::mouseMoveEvent(QMouseEvent* event) {
    QPointF worldPos = screenToWorld(event->pos());
    
    if (resizing_) {
        QPointF delta = worldPos - dragStart_;
        QSizeF newSize = originalSize_ + QSizeF(delta.x(), delta.y());
        
        // Ensure minimum size
        if (newSize.width() > 0.1 && newSize.height() > 0.1) {
            queryRegion_ = QRectF(originalTopLeft_, newSize);
            updateSearch();
            update();
        }
    }
    else if (dragging_) {
        QPointF delta = worldPos - dragStart_;
        queryRegion_.moveTopLeft(originalTopLeft_ + delta);
        
        // Keep within bounds
        if (queryRegion_.left() < 0) queryRegion_.moveLeft(0);
        if (queryRegion_.top() < 0) queryRegion_.moveTop(0);
        if (queryRegion_.right() > WORLD_SIZE) queryRegion_.moveRight(WORLD_SIZE);
        if (queryRegion_.bottom() > WORLD_SIZE) queryRegion_.moveBottom(WORLD_SIZE);
        
        updateSearch();
        update();
    }
}

void RangeSearchWidget::mouseReleaseEvent(QMouseEvent* event) {
    if (event->button() == Qt::LeftButton) {
        dragging_ = false;
        resizing_ = false;
    }
}

void RangeSearchWidget::generateRandomPoints() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.5, WORLD_SIZE - 0.5);
    
    points_.clear();
    
    for (int i = 0; i < pointCount_; ++i) {
        double x = dis(gen);
        double y = dis(gen);
        points_.push_back(Point2D(x, y));
    }
    
    // Rebuild KD-tree
    kdTree_ = std::make_unique<KDTree<2>>(points_);
    updateSearch();
    update();
}

void RangeSearchWidget::setPointCount(int count) {
    pointCount_ = count;
}

void RangeSearchWidget::toggleShowTree() {
    showTreePartitions_ = !showTreePartitions_;
    update();
}

void RangeSearchWidget::updateSearch() {
    searchResults_.clear();
    
    if (!kdTree_ || kdTree_->empty()) return;
    
    // Create CGAL range query
    Point2D minPoint(queryRegion_.left(), queryRegion_.top());
    Point2D maxPoint(queryRegion_.right(), queryRegion_.bottom());
    
    Kernel<2> kernel;
    Iso_box2D box = kernel.construct_iso_box_d_object()(minPoint, maxPoint);
    
    searchResults_ = kdTree_->range_search(box);
}

// Main Window Implementation
RangeSearchMainWindow::RangeSearchMainWindow(QWidget* parent)
    : QMainWindow(parent)
{
    setWindowTitle("Range Search Visualization");
    setMinimumSize(1000, 700);
    
    // Create central widget
    QWidget* centralWidget = new QWidget;
    setCentralWidget(centralWidget);
    
    // Create layout
    QVBoxLayout* mainLayout = new QVBoxLayout(centralWidget);
    
    // Create control panel
    QHBoxLayout* controlLayout = new QHBoxLayout;
    
    QPushButton* generateBtn = new QPushButton("Generate Random Points");
    showTreeCheckBox_ = new QCheckBox("Show KD-Tree Partitions");
    
    pointCountSpinBox_ = new QSpinBox;
    pointCountSpinBox_->setRange(1, 10000);
    pointCountSpinBox_->setValue(20);
    
    controlLayout->addWidget(new QLabel("Points:"));
    controlLayout->addWidget(pointCountSpinBox_);
    controlLayout->addWidget(generateBtn);
    controlLayout->addWidget(showTreeCheckBox_);
    controlLayout->addStretch();
    
    visualizationWidget_ = new RangeSearchWidget;
    
    mainLayout->addLayout(controlLayout);
    mainLayout->addWidget(visualizationWidget_);
    
    // Connect signals
    connect(generateBtn, &QPushButton::clicked, visualizationWidget_, &RangeSearchWidget::generateRandomPoints);
    connect(showTreeCheckBox_, &QCheckBox::toggled, visualizationWidget_, &RangeSearchWidget::toggleShowTree);
    connect(pointCountSpinBox_, QOverload<int>::of(&QSpinBox::valueChanged), 
            visualizationWidget_, &RangeSearchWidget::setPointCount);
}
