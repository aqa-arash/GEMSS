// gemss_gui.cpp
// Standalone Qt/VTK GUI for GEMSS multisphere-cpp
// Requirements: Qt5/6, VTK, Eigen, GEMSS headers and library

#include <QApplication>
#include <QMainWindow>
#include <QPushButton>
#include <QSlider>
#include <QLabel>
#include <QLineEdit>
#include <QFormLayout>
#include <QTextEdit>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QWidget>
#include <QGroupBox>

#include <QDoubleValidator>
#include <QSpinBox>
#include <QCheckBox>

#include <QVTKOpenGLNativeWidget.h>

#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSTLReader.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkSphereSource.h>
#include <vtkAppendPolyData.h>
#include <vtkSmartPointer.h>

#include <vtkNamedColors.h>
#include <vtkCamera.h>
#include <vtkAxesActor.h>
#include <vtkCubeSource.h>

#include "GEMSS/GEMSS-interface.h"

class GEMSSWindow : public QMainWindow {
    Q_OBJECT
public:
    GEMSSWindow(QWidget *parent = nullptr) : QMainWindow(parent), axes(nullptr), cubeActor(nullptr) {
        setWindowTitle("GEMSS Multisphere Visualizer");
        resize(1400, 900);
        QWidget *central = new QWidget(this);
        QHBoxLayout *mainLayout = new QHBoxLayout(central);
        QVBoxLayout *controlsLayout = new QVBoxLayout();
        // File load
        QPushButton *loadBtn = new QPushButton("Load STL");
        connect(loadBtn, &QPushButton::clicked, this, &GEMSSWindow::onLoadSTL);
        controlsLayout->addWidget(loadBtn);
        // Opacity
        QLabel *opacityLabel = new QLabel("Mesh Opacity");
        QSlider *opacitySlider = new QSlider(Qt::Horizontal);
        opacitySlider->setRange(1, 100);
        opacitySlider->setValue(100);
        connect(opacitySlider, &QSlider::valueChanged, this, &GEMSSWindow::onOpacityChanged);
        controlsLayout->addWidget(opacityLabel);
        controlsLayout->addWidget(opacitySlider);
        // Config
        QGroupBox *configBox = new QGroupBox("MultisphereConfig");
        QFormLayout *configForm = new QFormLayout();
        // --- Add all config parameters ---
        divSpin = new QSpinBox(); divSpin->setRange(1, 1000); divSpin->setValue(100);
        configForm->addRow("div", divSpin);
        paddingSpin = new QSpinBox(); paddingSpin->setRange(0, 100); paddingSpin->setValue(2);
        configForm->addRow("padding", paddingSpin);
        minRadiusRealEdit = new QLineEdit("0.0");
        minRadiusRealEdit->setValidator(new QDoubleValidator(0, 1e6, 4, this));
        configForm->addRow("minimum_radius_real", minRadiusRealEdit);
        confineMeshBox = new QCheckBox();
        configForm->addRow("confine_mesh", confineMeshBox);
        minCenterDistEdit = new QLineEdit("0.5");
        minCenterDistEdit->setValidator(new QDoubleValidator(0, 1e6, 4, this));
        configForm->addRow("min_center_distance_rel", minCenterDistEdit);
        minRadiusVoxSpin = new QSpinBox(); minRadiusVoxSpin->setRange(1, 100); minRadiusVoxSpin->setValue(2);
        configForm->addRow("min_radius_vox", minRadiusVoxSpin);
        precisionTargetEdit = new QLineEdit("1.0");
        precisionTargetEdit->setValidator(new QDoubleValidator(0, 1.0, 4, this));
        configForm->addRow("precision_target", precisionTargetEdit);
        maxSpheresSpin = new QSpinBox(); maxSpheresSpin->setRange(0, 1000000); maxSpheresSpin->setValue(0);
        configForm->addRow("max_spheres", maxSpheresSpin);
        computePhysicsSpin = new QSpinBox(); computePhysicsSpin->setRange(0, 2); computePhysicsSpin->setValue(2);
        configForm->addRow("compute_physics", computePhysicsSpin);
        pruneIsolatedBox = new QCheckBox();
        configForm->addRow("prune_isolated_spheres", pruneIsolatedBox);
        showProgressBox = new QCheckBox(); showProgressBox->setChecked(true);
        configForm->addRow("show_progress", showProgressBox);
        radiusOffsetEdit = new QLineEdit("0.87");
        radiusOffsetEdit->setValidator(new QDoubleValidator(0, 1e6, 4, this));
        configForm->addRow("radius_offset_vox", radiusOffsetEdit);
        searchWindowSpin = new QSpinBox(); searchWindowSpin->setRange(1, 100); searchWindowSpin->setValue(2);
        configForm->addRow("search_window", searchWindowSpin);
        persistenceSpin = new QSpinBox(); persistenceSpin->setRange(1, 100); persistenceSpin->setValue(2);
        configForm->addRow("persistence", persistenceSpin);
        configBox->setLayout(configForm);
        controlsLayout->addWidget(configBox);
        // Run button
        QPushButton *runBtn = new QPushButton("Run Multisphere");
        connect(runBtn, &QPushButton::clicked, this, &GEMSSWindow::onRunMultisphere);
        controlsLayout->addWidget(runBtn);
        // Info output
        infoBox = new QTextEdit(); infoBox->setReadOnly(true);
        controlsLayout->addWidget(new QLabel("Sphere Pack Info"));
        controlsLayout->addWidget(infoBox);
        controlsLayout->addStretch();
        // VTK
        vtkWidget = new QVTKOpenGLNativeWidget(this);
        renderWindow = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
        vtkWidget->setRenderWindow(renderWindow);
        renderer = vtkSmartPointer<vtkRenderer>::New();
        renderWindow->AddRenderer(renderer);

        // Add XYZ axes
        axes = vtkSmartPointer<vtkAxesActor>::New();
        axes->SetTotalLength(20.0, 20.0, 20.0); // Length of axes
        axes->SetShaftTypeToCylinder();
        axes->SetCylinderRadius(0.05);
        axes->SetAxisLabels(1);
        renderer->AddActor(axes);

        // Add a unit cube at the origin for scale reference
        vtkSmartPointer<vtkCubeSource> unitCube = vtkSmartPointer<vtkCubeSource>::New();
        unitCube->SetXLength(1.0);
        unitCube->SetYLength(1.0);
        unitCube->SetZLength(1.0);
        unitCube->SetCenter(0.5, 0.5, 0.5); // So it sits at origin
        vtkSmartPointer<vtkPolyDataMapper> cubeMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        cubeMapper->SetInputConnection(unitCube->GetOutputPort());
        cubeActor = vtkSmartPointer<vtkActor>::New();
        cubeActor->SetMapper(cubeMapper);
        cubeActor->GetProperty()->SetColor(0.2, 0.2, 0.2);
        cubeActor->GetProperty()->SetOpacity(0.3);
        renderer->AddActor(cubeActor);
        mainLayout->addLayout(controlsLayout, 1);
        mainLayout->addWidget(vtkWidget, 3);
        setCentralWidget(central);
    }

private slots:
    void onLoadSTL() {
        QString fname = QFileDialog::getOpenFileName(this, "Open STL File", "", "STL Files (*.stl)");
        if (fname.isEmpty()) return;
        currentSTL = fname;
        vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
        reader->SetFileName(fname.toStdString().c_str());
        reader->Update();
        meshActor = vtkSmartPointer<vtkActor>::New();
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(reader->GetOutputPort());
        meshActor->SetMapper(mapper);
        meshActor->GetProperty()->SetOpacity(1.0);
        renderer->RemoveAllViewProps();
        renderer->AddActor(meshActor);
        if (axes) renderer->AddActor(axes);
        if (cubeActor) renderer->AddActor(cubeActor);
        renderer->ResetCamera();
        renderWindow->Render();
    }
    void onOpacityChanged(int value) {
        if (meshActor) {
            meshActor->GetProperty()->SetOpacity(value / 100.0);
            renderWindow->Render();
        }
    }
    void onRunMultisphere() {
        if (currentSTL.isEmpty()) return;
        // --- Load mesh using GEMSS ---
        GEMSS::FastMesh mesh = GEMSS::load_mesh_fast(currentSTL.toStdString());
        GEMSS::MultisphereConfig config;
        config.div = divSpin->value();
        config.padding = paddingSpin->value();
        config.minimum_radius_real = minRadiusRealEdit->text().toFloat();
        config.confine_mesh = confineMeshBox->isChecked();
        config.min_center_distance_rel = minCenterDistEdit->text().toFloat();
        config.min_radius_vox = minRadiusVoxSpin->value();
        config.precision_target = precisionTargetEdit->text().toFloat();
        config.max_spheres = maxSpheresSpin->value();
        config.compute_physics = computePhysicsSpin->value();
        config.prune_isolated_spheres = pruneIsolatedBox->isChecked();
        config.show_progress = showProgressBox->isChecked();
        config.radius_offset_vox = radiusOffsetEdit->text().toFloat();
        config.search_window = searchWindowSpin->value();
        config.persistence = persistenceSpin->value();
        // --- Run reconstruction ---
        GEMSS::SpherePack pack = GEMSS::multisphere_from_mesh(mesh, config);
        // --- Visualize spheres ---
        // Color code spheres by radius
        float min_r = pack.radii.size() > 0 ? pack.radii.minCoeff() : 0.0f;
        float max_r = pack.radii.size() > 0 ? pack.radii.maxCoeff() : 1.0f;
        std::vector<vtkSmartPointer<vtkActor>> sphereActors;
        for (int i = 0; i < pack.centers.rows(); ++i) {
            vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
            sphere->SetCenter(pack.centers(i,0), pack.centers(i,1), pack.centers(i,2));
            sphere->SetRadius(pack.radii(i));
            sphere->SetThetaResolution(16);
            sphere->SetPhiResolution(16);
            sphere->Update();
            vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper->SetInputConnection(sphere->GetOutputPort());
            vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
            actor->SetMapper(mapper);
            // Normalize radius to [0,1] for colormap
            float norm = (pack.radii(i) - min_r) / (max_r - min_r + 1e-6f);
            // Simple blue-to-red colormap
            double r = norm;
            double g = 0.2;
            double b = 1.0 - norm;
            actor->GetProperty()->SetColor(r, g, b);
            actor->GetProperty()->SetOpacity(0.5);
            sphereActors.push_back(actor);
        }
        renderer->RemoveAllViewProps();
        renderer->AddActor(meshActor);
        for (auto &actor : sphereActors) {
            renderer->AddActor(actor);
        }
        if (axes) renderer->AddActor(axes);
        if (cubeActor) renderer->AddActor(cubeActor);
        renderer->ResetCamera();
        renderWindow->Render();
        // --- Info ---
        QString info = QString("Spheres: %1\nPrecision: %2\nVolume: %3\nBounding radius: %4\nCoM: [%5, %6, %7]\n")
            .arg(pack.num_spheres())
            .arg(pack.precision)
            .arg(pack.volume)
            .arg(pack.bounding_radius)
            .arg(pack.center_of_mass.x()).arg(pack.center_of_mass.y()).arg(pack.center_of_mass.z());

        // Moments of inertia (principal moments)
        info += QString("Principal moments: [%1, %2, %3]\n")
            .arg(pack.principal_moments(0))
            .arg(pack.principal_moments(1))
            .arg(pack.principal_moments(2));

        // Principal axes (as columns of the matrix)
        info += "Principal axes (columns):\n";
        for (int col = 0; col < 3; ++col) {
            info += QString("  [%1, %2, %3]\n")
                .arg(pack.principal_axes(0, col))
                .arg(pack.principal_axes(1, col))
                .arg(pack.principal_axes(2, col));
        }

        infoBox->setText(info);
    }
private:
    // Config widgets
    QSpinBox *divSpin, *paddingSpin, *minRadiusVoxSpin, *maxSpheresSpin, *computePhysicsSpin, *searchWindowSpin, *persistenceSpin;
    QLineEdit *minRadiusRealEdit, *minCenterDistEdit, *precisionTargetEdit, *radiusOffsetEdit;
    QCheckBox *confineMeshBox, *pruneIsolatedBox, *showProgressBox;
    QTextEdit *infoBox;
    // VTK
    QVTKOpenGLNativeWidget *vtkWidget;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow;
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkActor> meshActor;
    vtkSmartPointer<vtkAxesActor> axes;
    vtkSmartPointer<vtkActor> cubeActor;
    QString currentSTL;
};

#include <QVTKOpenGLNativeWidget.h>

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    GEMSSWindow win;
    win.show();
    return app.exec();
}

#include "gemss_gui.moc"
