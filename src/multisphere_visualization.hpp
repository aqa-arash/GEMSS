#ifndef MULTISPHERE_VISUALIZATION_HPP
#define MULTISPHERE_VISUALIZATION_HPP

#ifdef HAVE_VTK

#include "multisphere_datatypes.hpp"

#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkNamedColors.h>
#include <vtkLookupTable.h>
#include <vtkColorSeries.h>
#include <vtkAxesActor.h>
#include <vtkOrientationMarkerWidget.h>

void plot_sphere_pack(const SpherePack& sp, 
                      double opacity = 1.0, 
                      int phi_res = 50, 
                      int theta_res = 50,
                      std::string background = "white",
                      bool show_axes = true) {
    
    auto renderer = vtkSmartPointer<vtkRenderer>::New();
    auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    auto interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderWindow);

    auto colors = vtkSmartPointer<vtkNamedColors>::New();
    renderer->SetBackground(colors->GetColor3d(background).GetData());

    // Create a LookupTable for the colormap (Viridis style)
    auto lut = vtkSmartPointer<vtkLookupTable>::New();
    lut->SetNumberOfColors(256);
    lut->SetHueRange(0.667, 0.0); // Blue to Red
    lut->Build();

    for (int i = 0; i < sp.num_spheres(); ++i) {
        auto sphereSource = vtkSmartPointer<vtkSphereSource>::New();
        sphereSource->SetCenter(sp.centers(i, 0), sp.centers(i, 1), sp.centers(i, 2));
        sphereSource->SetRadius(sp.radii(i));
        sphereSource->SetPhiResolution(phi_res);
        sphereSource->SetThetaResolution(theta_res);

        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(sphereSource->GetOutputPort());

        auto actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetOpacity(opacity);
        
        // Color distribution t in [0, 1]
        double t = (double)i / std::max((size_t)1, sp.num_spheres() - 1);
        double rgb[3];
        lut->GetColor(t, rgb);
        actor->GetProperty()->SetColor(rgb);

        renderer->AddActor(actor);
    }

    if (show_axes) {
        auto axes = vtkSmartPointer<vtkAxesActor>::New();
        auto widget = vtkSmartPointer<vtkOrientationMarkerWidget>::New();
        widget->SetOrientationMarker(axes);
        widget->SetInteractor(interactor);
        widget->SetEnabled(1);
        widget->InteractiveOn();
    }

    renderWindow->SetWindowName("SpherePack Visualization");
    renderWindow->Render();
    interactor->Start();
}


void plot_mesh(const std::vector<Eigen::Vector3f>& vertices, 
               const std::vector<std::array<int, 3>>& faces,
               double opacity = 1.0,
               std::string color = "lightgray") {

    auto points = vtkSmartPointer<vtkPoints>::New();
    for (const auto& v : vertices) points->InsertNextPoint(v.x(), v.y(), v.z());

    auto vtk_faces = vtkSmartPointer<vtkCellArray>::New();
    for (const auto& f : faces) {
        vtkIdType cell[3] = {f[0], f[1], f[2]};
        vtk_faces->InsertNextCell(3, cell);
    }

    auto polyData = vtkSmartPointer<vtkPolyData>::New();
    polyData->SetPoints(points);
    polyData->SetPolys(vtk_faces);

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(polyData);

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    
    auto colors = vtkSmartPointer<vtkNamedColors>::New();
    actor->GetProperty()->SetColor(colors->GetColor3d(color).GetData());
    actor->GetProperty()->SetOpacity(opacity);
    actor->GetProperty()->SetInterpolationToPointLighting(); // Smooth shading

    auto renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->AddActor(actor);
    renderer->SetBackground(colors->GetColor3d("white").GetData());

    auto renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
    renderWindow->AddRenderer(renderer);
    auto interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(renderWindow);

    renderWindow->Render();
    interactor->Start();
}

#endif // for the vtk flag 
#endif // for header file 