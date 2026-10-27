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
#include <QProgressBar>
#include <QtConcurrent>
#include <QFutureWatcher>

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

#include <iostream>
#include <streambuf>
#include <exception>

#include "GEMSS/GEMSS-interface.h"

class StdRedirector : public QObject, public std::streambuf {
    Q_OBJECT
public:
    StdRedirector(std::ostream& stream) : stream(stream) {
        old_buf = stream.rdbuf(this);
    }
    ~StdRedirector() {
        stream.rdbuf(old_buf);
    }

signals:
    void textReceived(QString text);

protected:
    virtual int_type overflow(int_type v) override {
        if (v == '\n') {
            emit textReceived(QString::fromStdString(currentLine));
            currentLine.clear();
        } else {
            currentLine += static_cast<char>(v);
        }
        return v;
    }

private:
    std::ostream& stream;
    std::streambuf* old_buf;
    std::string currentLine;
};


class GEMSSWindow : public QMainWindow {
    Q_OBJECT
public:
    GEMSSWindow(QWidget *parent = nullptr) : QMainWindow(parent), axes(nullptr), cubeActor(nullptr) {
        setWindowTitle("GEMSS Multisphere Visualizer");
        resize(1400, 900);
        QWidget *central = new QWidget(this);
        QHBoxLayout *mainLayout = new QHBoxLayout(central);
        QVBoxLayout *controlsLayout = new QVBoxLayout();
        
        // File I/O layout
        QHBoxLayout *fileLayout = new QHBoxLayout();
        QPushButton *loadBtn = new QPushButton("Load STL");
        connect(loadBtn, &QPushButton::clicked, this, &GEMSSWindow::onLoadSTL);
        fileLayout->addWidget(loadBtn);

        importJsonBtn = new QPushButton("Import JSON");
        connect(importJsonBtn, &QPushButton::clicked, this, &GEMSSWindow::onImportJSON);
        fileLayout->addWidget(importJsonBtn);

        exportJsonBtn = new QPushButton("Export JSON");
        exportJsonBtn->setEnabled(false);
        connect(exportJsonBtn, &QPushButton::clicked, this, &GEMSSWindow::onExportJSON);
        fileLayout->addWidget(exportJsonBtn);

        controlsLayout->addLayout(fileLayout);

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
        runBtn = new QPushButton("Run Multisphere");
        connect(runBtn, &QPushButton::clicked, this, &GEMSSWindow::onRunMultisphere);        
        controlsLayout->addWidget(runBtn);
        
        // Initialize the Progress Bar
        progressBar = new QProgressBar();
        progressBar->setRange(0, 100);
        progressBar->setValue(0);
        progressBar->setTextVisible(false); // Keeps it clean
        controlsLayout->addWidget(progressBar);
        
        // Info output
        infoBox = new QTextEdit(); infoBox->setReadOnly(true);
        controlsLayout->addWidget(new QLabel("Sphere Pack Info"));
        controlsLayout->addWidget(infoBox);
        controlsLayout->addStretch();
        
        // Create the redirector for std::cout
        redirector = new StdRedirector(std::cout);
        
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
        
        // setup watcher
        watcher = new QFutureWatcher<GEMSS::SpherePack>(this);
        connect(watcher, &QFutureWatcherBase::finished, this, &GEMSSWindow::onProcessingFinished);
        
        // Connect the redirector to the text update slot
        connect(redirector, &StdRedirector::textReceived, 
            this, &GEMSSWindow::onTextReceived, Qt::QueuedConnection);
    }

private slots:
    void onTextReceived(QString text) {
        infoBox->append(text);
    }
    
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

    void onImportJSON() {
        QString fname = QFileDialog::getOpenFileName(this, "Import JSON File", "", "JSON Files (*.json)");
        if (fname.isEmpty()) return;
        
        try {
            infoBox->clear();
            infoBox->append("<b>--- Loading SpherePack from JSON ---</b>");
            currentPack = GEMSS::import_sphere_pack(fname.toStdString());
            exportJsonBtn->setEnabled(true);
            displaySpherePack();
        } catch (const std::exception& e) {
            infoBox->append(QString("<font color='red'>Failed to load JSON: %1</font>").arg(e.what()));
        }
    }

    void onExportJSON() {
        if (currentPack.num_spheres() == 0) return;
        QString fname = QFileDialog::getSaveFileName(this, "Export JSON File", "", "JSON Files (*.json)");
        if (fname.isEmpty()) return;
        
        try {
            GEMSS::export_sphere_pack_json(currentPack, fname.toStdString());
            infoBox->append(QString("<b>Exported to:</b> %1").arg(fname));
        } catch (const std::exception& e) {
            infoBox->append(QString("<font color='red'>Failed to export JSON: %1</font>").arg(e.what()));
        }
    }

    void onOpacityChanged(int value) {
        if (meshActor) {
            meshActor->GetProperty()->SetOpacity(value / 100.0);
            renderWindow->Render();
        }
    }

    void onRunMultisphere() {
        if (currentSTL.isEmpty()) return;
        infoBox->clear();
        infoBox->append("<b>--- Starting Multisphere Calculation ---</b>");
        runBtn->setEnabled(false);
        exportJsonBtn->setEnabled(false);
        progressBar->setRange(0, 0);
        infoBox->setText("Processing... please wait.");

        GEMSS::STLMesh mesh = GEMSS::load_mesh(currentSTL.toStdString());
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
        
        QFuture<GEMSS::SpherePack> future = QtConcurrent::run(GEMSS::multisphere_from_mesh, mesh, config);
        watcher->setFuture(future);        
    }

    void onProcessingFinished() { 
        currentPack = watcher->result();
        progressBar->setRange(0, 100);
        progressBar->setValue(100);
        runBtn->setEnabled(true);
        exportJsonBtn->setEnabled(true);
        
        displaySpherePack();
    }

    void displaySpherePack() {
        float min_r = currentPack.radii.size() > 0 ? currentPack.radii.minCoeff() : 0.0f;
        float max_r = currentPack.radii.size() > 0 ? currentPack.radii.maxCoeff() : 1.0f;
        
        std::vector<vtkSmartPointer<vtkActor>> sphereActors;
        for (int i = 0; i < currentPack.centers.rows(); ++i) {
            vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
            sphere->SetCenter(currentPack.centers(i,0), currentPack.centers(i,1), currentPack.centers(i,2));
            sphere->SetRadius(currentPack.radii(i));
            sphere->SetThetaResolution(16);
            sphere->SetPhiResolution(16);
            sphere->Update();
            
            vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper->SetInputConnection(sphere->GetOutputPort());
            vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
            actor->SetMapper(mapper);
            
            float norm = (currentPack.radii(i) - min_r) / (max_r - min_r + 1e-6f);
            double r = norm;
            double g = 0.2;
            double b = 1.0 - norm;
            actor->GetProperty()->SetColor(r, g, b);
            actor->GetProperty()->SetOpacity(0.5);
            sphereActors.push_back(actor);
        }
        
        renderer->RemoveAllViewProps();
        if (meshActor) renderer->AddActor(meshActor);
        for (auto &actor : sphereActors) {
            renderer->AddActor(actor);
        }
        if (axes) renderer->AddActor(axes);
        if (cubeActor) renderer->AddActor(cubeActor);
        renderer->ResetCamera();
        renderWindow->Render();
        
        QString info = QString("Spheres: %1\nPrecision: %2\nMass: %3\nBounding radius: %4\nCoM: [%5, %6, %7]\n")
            .arg(currentPack.num_spheres())
            .arg(currentPack.precision)
            .arg(currentPack.mass)
            .arg(currentPack.bounding_radius)
            .arg(currentPack.center_of_mass.x()).arg(currentPack.center_of_mass.y()).arg(currentPack.center_of_mass.z());

        info += QString("Principal moments: [%1, %2, %3]\n")
            .arg(currentPack.principal_moments(0))
            .arg(currentPack.principal_moments(1))
            .arg(currentPack.principal_moments(2));

        info += "Principal axes (columns):\n";
        for (int col = 0; col < 3; ++col) {
            info += QString("  [%1, %2, %3]\n")
                .arg(currentPack.principal_axes(0, col))
                .arg(currentPack.principal_axes(1, col))
                .arg(currentPack.principal_axes(2, col));
        }

        info += QString("Principal extents (L, I, S): [%1, %2, %3]\n")
            .arg(currentPack.principal_extents(0))
            .arg(currentPack.principal_extents(1))
            .arg(currentPack.principal_extents(2));

        infoBox->setText(info);
    }

private:
    QSpinBox *divSpin, *paddingSpin, *minRadiusVoxSpin, *maxSpheresSpin, *computePhysicsSpin, *searchWindowSpin, *persistenceSpin;
    QLineEdit *minRadiusRealEdit, *minCenterDistEdit, *precisionTargetEdit, *radiusOffsetEdit;
    QCheckBox *confineMeshBox, *pruneIsolatedBox, *showProgressBox;
    QTextEdit *infoBox;
    QProgressBar *progressBar;
    QPushButton *runBtn, *importJsonBtn, *exportJsonBtn;
    QFutureWatcher<GEMSS::SpherePack> *watcher;
    
    QVTKOpenGLNativeWidget *vtkWidget;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow;
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkActor> meshActor;
    vtkSmartPointer<vtkAxesActor> axes;
    vtkSmartPointer<vtkActor> cubeActor;
    
    QString currentSTL;
    GEMSS::SpherePack currentPack;
    StdRedirector* redirector;
};

#include <QVTKOpenGLNativeWidget.h>

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    GEMSSWindow win;
    win.show();
    return app.exec();
}

#include "gemss_gui.moc"
