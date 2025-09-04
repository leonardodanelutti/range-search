#pragma once

#include <QtWidgets>
#include <QApplication>
#include <QMainWindow>
#include <QWidget>
#include <QPainter>
#include <QMouseEvent>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QLabel>
#include <QSpinBox>
#include <QCheckBox>
#include <vector>
#include <random>
#include "geometry.h"
#include "algorithms.h"

class RangeSearchWidget : public QWidget {
    Q_OBJECT

public:
    RangeSearchWidget(QWidget* parent = nullptr);

public slots:
    void generateRandomPoints();
    void setPointCount(int count);
    void toggleShowTree();
    
protected:
    void paintEvent(QPaintEvent* event) override;
    void mousePressEvent(QMouseEvent* event) override;
    void mouseMoveEvent(QMouseEvent* event) override;
    void mouseReleaseEvent(QMouseEvent* event) override;

private:
    // Convert between screen coordinates and world coordinates
    QPointF screenToWorld(const QPoint& screen) const;
    QPoint worldToScreen(const QPointF& world) const;
    
    // Drawing functions
    void drawPoints(QPainter& painter);
    void drawQueryRegion(QPainter& painter);
    void drawKDTreePartitions(QPainter& painter);
    void drawResultPoints(QPainter& painter);
    
    // Helper function to draw KD-tree partitions recursively
    void drawTreePartitions(QPainter& painter, KDTree<2>::Node* node, 
                          double minX, double minY, double maxX, double maxY, int depth) const;
    
    // Update the range search results
    void updateSearch();
    
    // Data members
    std::vector<Point2D> points_;
    std::unique_ptr<KDTree<2>> kdTree_;
    std::vector<Point2D> searchResults_;
    
    // Query region (in world coordinates)
    QRectF queryRegion_;
    
    // Interaction state
    bool dragging_;
    bool resizing_;
    QPointF dragStart_;
    QPointF originalTopLeft_;
    QSizeF originalSize_;
    
    // Visualization options
    bool showTreePartitions_;
    int pointCount_;
    
    // Transform parameters
    double scale_;
    QPointF offset_;
    
    // Widget dimensions
    static const int MARGIN = 50;
    static constexpr double WORLD_SIZE = 10.0;
};

class RangeSearchMainWindow : public QMainWindow {
    Q_OBJECT

public:
    RangeSearchMainWindow(QWidget* parent = nullptr);

private:
    RangeSearchWidget* visualizationWidget_;
    QSpinBox* pointCountSpinBox_;
    QCheckBox* showTreeCheckBox_;
};
