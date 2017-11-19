#include "VTKViewer.h"
#include <QString>

//class MainWindow;

// Catch mouse and keyboard events
class MouseInteractorStyle : public vtkInteractorStyleTrackballCamera
{

public:
	static MouseInteractorStyle* New();
	vtkTypeMacro(MouseInteractorStyle, vtkInteractorStyleTrackballCamera);

	MouseInteractorStyle()
	{
		selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
		selectedActor = vtkSmartPointer<vtkActor>::New();
		has_selected_one_landmark = false;
		//count_visible_points = false;
	}

	virtual void OnLeftButtonDown()
	{
		vtkRenderWindowInteractor *rwi = this->Interactor;
		if (rwi->GetKeySym())
		{
			if (rwi->GetShiftKey())
			{
				// Get the location of the click (in window coordinates)
				int* pos = this->GetInteractor()->GetEventPosition();

				picker = vtkSmartPointer<vtkCellPicker>::New();
				picker->SetTolerance(0.0005);

				// Pick from this location.
				picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

				double* worldPosition = picker->GetPickPosition();
				//std::cout << "Cell id is: " << picker->GetCellId() << std::endl;

				if (picker->GetCellId() != -1)
				{
					if (!has_selected_one_landmark)
					{
						// pick one landmark
						ids = vtkSmartPointer<vtkIdTypeArray>::New();
						vtkIdType nearest_landmark_id = this_viewer->findNearestLandmark(picker->GetCellId());
						if (nearest_landmark_id == -1) // 3D out of distance
						{
							refresh();
							return;
						}
						ids->SetNumberOfComponents(1);
						ids->InsertNextValue(nearest_landmark_id);
						has_selected_one_landmark = true;
						this_viewer->enableFineAnnotationButton();
					}
					else
					{
						vtkIdType new_position = picker->GetCellId();
						has_selected_one_landmark = false;
						// set it to the right place
						ids = vtkSmartPointer<vtkIdTypeArray>::New();
						ids->SetNumberOfComponents(1);
						//ids->InsertNextValue(new_position);
						// change the point in VTKViewer class
						this_viewer->changePositionOfLastSelectedLandmark(new_position); 
						this_viewer->to2DButtonPressed();
//						this_viewer->disableCancelButton();
					}
					refresh();
				}
			}
		}

		// Forward events
		vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
	}

	virtual void OnLeftButtonUp()
	{
		vtkInteractorStyleTrackballCamera::OnLeftButtonUp();
		/*
		// count visible points
		if (count_visible_points)
		{
			this_viewer->countPoints();
		}
		*/
	}

	virtual void OnRightButtonDown()
	{
		;
	}

	virtual void OnKeyPress()
	{
		// Get the keypress
		vtkRenderWindowInteractor *rwi = this->Interactor;
		std::string key = rwi->GetKeySym();

		// Output the key that was pressed
		//qDebug() << QString::fromStdString(key);

		// Handle a "normal" key
		if (key.compare("c") == 0)
		{
			//std::cout << "The a key was pressed." << std::endl;
			this_viewer->toggleVisibility();
		}
		/*else if (key.compare("s") == 0)
		{
			this_viewer->countPoints();
		}*/
		else if (key.compare("Escape") == 0)
		{
			this_viewer->onCancel();
		}

		// Forward events
		vtkInteractorStyleTrackballCamera::OnKeyPress();
	}

	void setViewerCallback(VTKViewer* viewer)
	{
		this_viewer = viewer;
	}

	vtkSmartPointer<vtkPolyData> Data;

	void cancel()
	{
		if (has_selected_one_landmark)
		{
			has_selected_one_landmark = false;
			ids = vtkSmartPointer<vtkIdTypeArray>::New();
			refresh();
		}
	}

	// for getting visible points 
	//void SetVisibleFilter(vtkSmartPointer<vtkSelectVisiblePoints> vis) { this->visible_filter = vis; }

protected:
	// color the selected node
	void refresh()
	{
		selectionNode = vtkSmartPointer<vtkSelectionNode>::New();
		selectionNode->SetFieldType(vtkSelectionNode::CELL);
		selectionNode->SetContentType(vtkSelectionNode::INDICES);
		selectionNode->SetSelectionList(ids);

		selection = vtkSmartPointer<vtkSelection>::New();
		selection->AddNode(selectionNode);

		extractSelection = vtkSmartPointer<vtkExtractSelection>::New();

		extractSelection->SetInputData(0, this->Data);
		extractSelection->SetInputData(1, selection);
		extractSelection->Update();

		// In selection
		selected = vtkSmartPointer<vtkUnstructuredGrid>::New();
		selected->ShallowCopy(extractSelection->GetOutput());

		// color the selected cell:
		selectedMapper->SetInputData(selected);

		selectedActor->SetMapper(selectedMapper);
		selectedActor->GetProperty()->EdgeVisibilityOn();
		selectedActor->GetProperty()->SetEdgeColor(1, 1, 0); // selected cell goes green
		selectedActor->GetProperty()->SetLineWidth(5);

		this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selectedActor);		
		this_viewer->refreshVisibility();
	}
	
	vtkSmartPointer<vtkDataSetMapper> selectedMapper;
	vtkSmartPointer<vtkActor> selectedActor;
	vtkSmartPointer<vtkIdTypeArray> ids;
	vtkSmartPointer<vtkCellPicker> picker;
	vtkSmartPointer<vtkSelectionNode> selectionNode;
	vtkSmartPointer<vtkSelection> selection;
	vtkSmartPointer<vtkExtractSelection> extractSelection;
	vtkSmartPointer<vtkUnstructuredGrid> selected;
	VTKViewer* this_viewer;

	bool has_selected_one_landmark;

	// for testing the visibility of landmarks
	//vtkSmartPointer<vtkSelectVisiblePoints> visible_filter;
	//bool count_visible_points;
};

vtkStandardNewMacro(MouseInteractorStyle);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


VTKViewer::VTKViewer(QWidget *parent) :
	QWidget(parent)
{

	//parent_pointer = (QMainWindow*)parentWidget();

	has_chosen_one_annotation = false;
	show_mesh = true;
	show_popup = true;
	// if initialize here, cannot be always on top
	//hint_dialog = new hintDialog();
	now_chosen_cell_idx = -1;
	last_changed_landmark_idx = -1;
	now_chosen_id_in_array = -1;

	vertical_layout = new QVBoxLayout(this);
	horizontal_layout = new QHBoxLayout(this);
	button_layout_left = new QVBoxLayout(this);
	button_layout_right = new QVBoxLayout(this);

	qvtkWidget = new QVTKWidget(this);

	// Undo Button
	//undo_button = new QPushButton("Undo", this);
	//connect(undo_button, &QPushButton::clicked, this, &VTKViewer::debugSetMesh);
		
	correct_this_point_button = new QPushButton("Correct This Point", this);
	correct_this_point_button->setDisabled(true);
	connect(correct_this_point_button, &QPushButton::clicked, this, &VTKViewer::correctThisPointButtonPressed_);

	to_2d_button = new QPushButton(" 3D >> 2D ", this);
	connect(to_2d_button, &QPushButton::clicked, this, &VTKViewer::to2DButtonPressed_);
	to_2d_button->hide();

	//connect(cancel_button, &QPushButton::clicked, this, &VTKViewer::onCancel);
	//cancel_button->setDisabled(true);
	//cancel_button->hide();

	// Toggle mesh display
	switch_num_button = new QPushButton("Toggle Mesh Display", this);
	connect(switch_num_button, &QPushButton::clicked, this, &VTKViewer::toggleVisibility);
	switch_num_button->setDisabled(true);
	switch_num_button->hide();

	// Toggle floating window
	//switch_instruction_button = new QPushButton("Toggle Instructions", this);
	//connect(switch_instruction_button, &QPushButton::clicked, this, &VTKViewer::popUpInstruction);
	//switch_instruction_button->hide();

	// Propagate back to 2D
	//propagate_back_to_2D_button = new QPushButton("Propagate Back", this);
	//connect(propagate_back_to_2D_button, &QPushButton::clicked, this, &VTKViewer::popUpInstruction);

	status_bar = new QStatusBar(this);
	status_bar->setFixedHeight(20);
	status_bar->showMessage("Hints: \"Mouse Left Button\" : Move. \"Shift + Mouse Left Button\" : Select 3D Point.  \"C\" : Toggle Display Mode. \"Esc\" : Cancel Selection.");

	// Set up the QVTK part ///////////////////////////////////////////////////////////////////////////////////////////

	helper_mesh.reset(new pcl::PolygonMesh);

	poly_data = vtkSmartPointer<vtkPolyData>::New();
	pcl::VTKUtils::convertToVTK(*helper_mesh, poly_data);

	// Visualize
	mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(poly_data);
	mapper->Update();

	actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	renderer = vtkSmartPointer<vtkRenderer>::New();
	renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(qvtkWidget->GetRenderWindow());
	renderWindowInteractor->Initialize();
	renderer->AddActor(actor);

	// Set the custom stype to use for interaction.
	style = vtkSmartPointer<MouseInteractorStyle>::New();
	style->SetDefaultRenderer(renderer);
	style->Data = poly_data;
	style->setViewerCallback(this);

	/*
	// set up the visible-point counter
	visible_filter = vtkSmartPointer<vtkSelectVisiblePoints>::New();
	//visible_filter->SetInputConnection();
	//visible_filter->SetInputConnection(poly_data->)
	//visible_filter->SetInputData(poly_data);
	visible_filter->SetInputData(0, poly_data);
	//visible_filter->SetInputDataObject(poly_data);
	//visible_filter->SetInputData(poly_data);
	visible_filter->SetRenderer(renderer);
	visible_filter->SelectionWindowOff();
	visible_filter->Update();
	*/

	renderWindowInteractor->SetInteractorStyle(style);
	//style->SetVisibleFilter(visible_filter);
	
	renderer->SetBackground(0, 0, 0); // Black background

	renderer->ResetCamera();
	qvtkWidget->GetRenderWindow()->AddRenderer(renderer);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	qvtkWidget->update();

	horizontal_layout->addWidget(correct_this_point_button);
	horizontal_layout->addStretch();
	horizontal_layout->addWidget(to_2d_button);

	vertical_layout->addLayout(horizontal_layout);
	vertical_layout->addWidget(qvtkWidget);
	vertical_layout->addWidget(status_bar);
	//horizontal_layout->addWidget(switch_num_button);
	//horizontal_layout->addWidget(switch_instruction_button);

	vertical_layout->setAlignment(Qt::AlignTop);
	vertical_layout->setContentsMargins(0, 0, 0, 0);

	// for selection highlighting
	selection_node = vtkSmartPointer<vtkSelectionNode>::New();
	selection = vtkSmartPointer<vtkSelection>::New();
	extract_selection = vtkSmartPointer<vtkExtractSelection>::New();
	selected = vtkSmartPointer<vtkUnstructuredGrid>::New();
	selected_mapper = vtkSmartPointer<vtkDataSetMapper>::New();
	selected_actor = vtkSmartPointer<vtkActor>::New();
}

void VTKViewer::setMesh(std::string mesh_path)
{
	// cache
	now_mesh_path = mesh_path;
	//
	helper_mesh.reset(new pcl::PolygonMesh);
	pcl::io::load(mesh_path, (*helper_mesh));
	pcl::VTKUtils::convertToVTK((*helper_mesh), poly_data);
	mapper->SetInputData(poly_data);
	mapper->Update();
	style->Data = poly_data;
	renderer->ResetCamera();
	qvtkWidget->update();
	//
	correct_this_point_button->setDisabled(true);
	//
	emit meshPointerSignal(poly_data);
}

void VTKViewer::debugSetMesh()
{
	setMesh("C:\\Users\\shangxuanu\\Desktop\\image_original\\mesh.ply");
}

void VTKViewer::pickPointCallback(vtkIdType idx)
{
	now_chosen_cell_idx = idx;
}

vtkIdType VTKViewer::getNearestCell(float x, float y, float z)
{
	vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
	cellLocator->SetDataSet(poly_data);
	cellLocator->BuildLocator();

	double testPoint[3] = { x, y, z };

	//Find the closest points to TestPoint
	double closestPoint[3];//the coordinates of the closest point will be returned here
	double closestPointDist2; //the squared distance to the closest point will be returned here
	vtkIdType cell_id; //the cell id of the cell containing the closest point will be returned here
	int sub_id; //this is rarely used (in triangle strips only, I believe)
	cellLocator->FindClosestPoint(testPoint, closestPoint, cell_id, sub_id, closestPointDist2);

	return cell_id;
}

void VTKViewer::highlightAllLandmarks()
{
	selection_node = vtkSmartPointer<vtkSelectionNode>::New();
	selection_node->SetFieldType(vtkSelectionNode::CELL);
	selection_node->SetContentType(vtkSelectionNode::INDICES);
	selection_node->SetSelectionList(all_landmarks_cell_id);

	selection = vtkSmartPointer<vtkSelection>::New();
	selection->AddNode(selection_node);

	extract_selection = vtkSmartPointer<vtkExtractSelection>::New();

	extract_selection->SetInputData(0, poly_data);
	extract_selection->SetInputData(1, selection);
	extract_selection->Update();

	// In selection
	selected = vtkSmartPointer<vtkUnstructuredGrid>::New();
	selected->ShallowCopy(extract_selection->GetOutput());

	selected_mapper->SetInputData(selected);

	selected_actor->SetMapper(selected_mapper);
	selected_actor->GetProperty()->EdgeVisibilityOn();
	selected_actor->GetProperty()->SetEdgeColor(0.4, 1, 1); // cyan, this is the color for all nodes
	selected_actor->GetProperty()->SetLineWidth(3);
	
	renderer->AddActor(selected_actor);

	// add text
	highlightAllLandmarks_text();

	// add line
	highlightAllLandmarks_line();
	
	// set mesh & points to be shown
	//setInitialVisibility();
}

/*
void VTKViewer::setInitialVisibility()
{
	vtkPropCollection* props = renderer->GetViewProps(); //iterate through and set each visibility to 0
	props->InitTraversal();
	props->GetNextProp()->VisibilityOn();
	props->GetNextProp()->VisibilityOn();
	for (int i = 2; i <= 70; i++)
	{
		props->GetNextProp()->VisibilityOff();
	}
	//props->GetNextProp()->VisibilityOff();
	qvtkWidget->update();
}
*/

void VTKViewer::toggleVisibility()
{
	/*if (!switch_num_button->isEnabled())
	{
		return;
	}*/

	if (show_mesh)
	{
		show_mesh = false;
	}
	else
	{
		show_mesh = true;
	}
	refreshVisibility();
}

void VTKViewer::refreshVisibility()
{
	// helper function
	//printVisibility();

	//
	vtkPropCollection* props = renderer->GetViewProps(); //iterate through and set each visibility to what we want
	props->InitTraversal();

	if (!show_mesh) // show line & text
	{
		props->GetNextProp()->VisibilityOff(); // mesh
		props->GetNextProp()->VisibilityOff(); // landmark
		for (int i = 2; i <= 2+ TOTAL_POINTS_IN_SINGLE_VIEW; i++)
		{
			props->GetNextProp()->VisibilityOn(); // 68 texts + 1 line
		}
		if (props->GetNumberOfItems() > 3 + TOTAL_POINTS_IN_SINGLE_VIEW)
		{
			props->GetNextProp()->VisibilityOff();// selected point
		}
	}
	else // show mesh & landmark
	{
		qDebug() << props->GetNumberOfItems();

		props->GetNextProp()->VisibilityOn();
		props->GetNextProp()->VisibilityOn();
		for (int i = 2; i <= 2 + TOTAL_POINTS_IN_SINGLE_VIEW; i++)
		{
			props->GetNextProp()->VisibilityOff();
		}
		if (props->GetNumberOfItems() > 3 + TOTAL_POINTS_IN_SINGLE_VIEW)
		{
			props->GetNextProp()->VisibilityOn();// selected point
		}
	}
	qvtkWidget->update();
}

// this is a helper function
void VTKViewer::printVisibility()
{
	vtkPropCollection* props = renderer->GetViewProps(); //iterate through and set each visibility to what we want
	props->InitTraversal();

	for (int i = 0; i < props->GetNumberOfItems(); i++)
	{
		qDebug() << i;
		std::string a = props->GetNextProp()->GetClassName();
		qDebug() << QString::fromStdString(a);
	}
}

void VTKViewer::highlightAllLandmarks_text()
{
	for (int i = 0; i < TOTAL_POINTS_IN_SINGLE_VIEW; i++)
	{
		// Create some text
		text_source = vtkSmartPointer<vtkVectorText>::New();
		text_source->SetText(std::to_string(i+1).c_str());

		//double* bounds = text_source->GetOutput()->GetBounds();
		//transform the polydata to be centered over the pick position
		//double center[3] = { 0.5*(bounds[1] + bounds[0]), 0.5*(bounds[3] + bounds[2]), 0.0 };

		double p[3];
		getCenterPointFromCellId(all_landmarks_cell_id->GetValue(i), p);

		// Create a mapper
		text_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		text_mapper->SetInputConnection(text_source->GetOutputPort());
	
		// Create a subclass of vtkActor: a vtkFollower that remains facing the camera
		follower = vtkSmartPointer<vtkFollower>::New();
		follower->SetMapper(text_mapper);
		follower->GetProperty()->SetColor(0.4, 1, 1); // cyan
		follower->SetPosition(p);
		follower->SetCamera(renderer->GetActiveCamera());

		// Add the actor to the scene
		renderer->AddActor(follower);
		followers.push_back(follower);
	}
}

void VTKViewer::highlightAllLandmarks_line()
{
	// Create the polydata where we will store all the geometric data
	lines_poly_data = vtkSmartPointer<vtkPolyData>::New();

	// Create a vtkPoints container and store the points in it
	pts = vtkSmartPointer<vtkPoints>::New();
	for (int i = 0; i < TOTAL_POINTS_IN_SINGLE_VIEW; i++)
	{
		double p[3];
		// get point from poly_data
		getCenterPointFromCellId(all_landmarks_cell_id->GetValue(i), p);
		//
		pts->InsertNextPoint(p);
	}

	// Add the points to the polydata container
	lines_poly_data->SetPoints(pts);
	lines = vtkSmartPointer<vtkCellArray>::New();
	
	for (auto it = connected_points.begin(); it != connected_points.end(); it++)
	{
		// Create the first line (between Origin and P0)
		vtkSmartPointer<vtkLine> line0 = vtkSmartPointer<vtkLine>::New();
		line0->GetPointIds()->SetId(0, it->first - 1); // the second 0 is the index of the Origin in linesPolyData's points
		line0->GetPointIds()->SetId(1, it->second - 1); // the second 1 is the index of P0 in linesPolyData's points
		
		lines->InsertNextCell(line0);
	}
	
		// Add the lines to the polydata container
	lines_poly_data->SetLines(lines);
	
	// Create two colors - one for each line
	unsigned char green[3] = { 0, 255, 0 };
	
	// Create a vtkUnsignedCharArray container and store the colors in it
	colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	for (int i = 0; i < connected_points.size(); i++)
	{
		colors->InsertNextTupleValue(green);
	}

	// Color the lines.
	// SetScalars() automatically associates the values in the data array passed as parameter
	// to the elements in the same indices of the cell data array on which it is called.
	// This means the first component (red) of the colors array
	// is matched with the first component of the cell array (line 0)
	// and the second component (green) of the colors array
	// is matched with the second component of the cell array (line 1)
	lines_poly_data->GetCellData()->SetScalars(colors);

	// Setup the visualization pipeline
	line_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	
	line_mapper->SetInputData(lines_poly_data);

	line_actor = vtkSmartPointer<vtkActor>::New();
	line_actor->SetMapper(line_mapper);

	renderer->AddActor(line_actor);
}

std::vector<pts_3d_conf> VTKViewer::setAndShowLandmarks(std::vector<pts_3d_conf> points)
{
	// should clear all landmarks first
	removeRedundantActors();
	//printVisibility();

	// then refresh all the landmarks
	std::vector<pts_3d_conf> result = setLandmarks(points);
	highlightAllLandmarks();
	refreshVisibility();
	switch_num_button->setEnabled(true);
	//initializePopUpInstruction();
	return result;
}

void VTKViewer::removeRedundantActors()
{
	vtkPropCollection* props = renderer->GetViewProps(); //iterate through and set each visibility to what we want
	props->InitTraversal();

	while (props->GetNumberOfItems() > 2)
	{
		renderer->RemoveActor(props->GetLastProp());
	}
	/*
	for (int i = 0; i < props->GetNumberOfItems(); i++)
	{
		if (i >=2)
		{
			renderer->RemoveActor(props->GetNextProp());
		}
	}
	*/
}

std::vector<pts_3d_conf> VTKViewer::setLandmarks(std::vector<pts_3d_conf> points)
{
	all_landmarks_cell_id = vtkSmartPointer<vtkIdTypeArray>::New();
	std::vector<pts_3d_conf> result;
	for (int i = 0; i < TOTAL_POINTS_IN_SINGLE_VIEW; i++)
	{
		// search for the nearest cell by function
		vtkIdType now_id = getNearestCell(points[i].x, points[i].y, points[i].z);
		all_landmarks_cell_id->InsertNextValue(now_id);
		double p[3];
		poly_data->GetCell(now_id)->GetPoints()->GetPoint(0,p);
		pts_3d_conf this_point = pts_3d_conf(p[0], p[1], p[2]);
		result.push_back(this_point);
	}
	return result;
}

/*
void VTKViewer::initializePopUpInstruction()
{
	hint_dialog = new hintDialog();
	hint_dialog->show();
}
*/

/*
void VTKViewer::popUpInstruction()
{
	if (!hint_dialog)
	{
		return;
	}
	if (!show_popup)
	{
		hint_dialog->show();
		show_popup = true;
	}
	else
	{
		hint_dialog->hide();
		show_popup = false;
	}
}
*/

/*
// used by mainwineow::close()
void VTKViewer::closeDialog()
{
	if (hint_dialog)
	{
		hint_dialog->close();
	}
}
*/

vtkIdType VTKViewer::findNearestLandmark(vtkIdType now_idx)
{
	std::vector<double> distances;
	for (int i = 0; i < TOTAL_POINTS_IN_SINGLE_VIEW; i++)
	{
		// count distance
		double distance = 0.0;
		
		for (int j = 0; j <= 2; j++)
		{
			vtkPoints* pts1 = poly_data->GetCell(now_idx)->GetPoints();
			double cords1[3];
			pts1->GetPoint(j, cords1);
			
			vtkIdType another_id = all_landmarks_cell_id->GetValue(i);
			vtkPoints* pts2 = poly_data->GetCell(another_id)->GetPoints();
			double cords2[3];
			pts2->GetPoint(j, cords2);
			
			distance += sqrt(vtkMath::Distance2BetweenPoints(cords1, cords2));
		}
		distances.push_back(distance);
	}
	double smallest_distance = INT_MAX;
	int smallest_distance_index = -1;
	for (int i = 0; i < TOTAL_POINTS_IN_SINGLE_VIEW; i++)
	{
		if (distances[i] < smallest_distance)
		{
			smallest_distance = distances[i];
			smallest_distance_index = i;
		}
	}
	if (smallest_distance < MAX_SELECT_DISTANCE_3D)
	{
		now_chosen_id_in_array = smallest_distance_index;
		//hint_dialog->setHint(smallest_distance_index);
		emit setHint(smallest_distance_index);
		emit setHint3D(smallest_distance_index);
		return all_landmarks_cell_id->GetValue(smallest_distance_index);
	}
	else
	{
		return -1;
	}
}

void VTKViewer::onCancel()
{
	style->cancel();
	qvtkWidget->update();
	//cancel_button->setDisabled(true);
	//hint_dialog->cancel();
	correct_this_point_button->setDisabled(true);
	emit cancelHint();
	emit cancelHint3D();
}

// this is actually the refresh function of this class
void VTKViewer::changePositionOfLastSelectedLandmark(vtkIdType new_cell_id)
{
	// change position in original cell array
	all_landmarks_cell_id->SetValue(now_chosen_id_in_array, new_cell_id);
	// update landmark shown
	selection_node->SetSelectionList(all_landmarks_cell_id);
	extract_selection->Update();
	selected->DeepCopy(extract_selection->GetOutput());
	selected_mapper->Update();
	// update text /////////////////////////////////////////////////////////////////////////////////
	double p[3];
	getCenterPointFromCellId(all_landmarks_cell_id->GetValue(now_chosen_id_in_array), p);
	followers[now_chosen_id_in_array]->SetPosition(p);
	// update line /////////////////////////////////////////////////////////////////////////////////
	// Create a vtkPoints container and store the points in it
	pts = vtkSmartPointer<vtkPoints>::New();
	for (int i = 0; i < TOTAL_POINTS_IN_SINGLE_VIEW; i++)
	{
		double p[3];
		// get point from poly_data
		getCenterPointFromCellId(all_landmarks_cell_id->GetValue(i), p);
		//
		pts->InsertNextPoint(p);
	}

	// Add the points to the polydata container
	lines_poly_data->SetPoints(pts);
	lines = vtkSmartPointer<vtkCellArray>::New();

	for (auto it = connected_points.begin(); it != connected_points.end(); it++)
	{
		// Create the first line (between Origin and P0)
		vtkSmartPointer<vtkLine> line0 = vtkSmartPointer<vtkLine>::New();
		line0->GetPointIds()->SetId(0, it->first - 1); // the second 0 is the index of the Origin in linesPolyData's points
		line0->GetPointIds()->SetId(1, it->second - 1); // the second 1 is the index of P0 in linesPolyData's points

		lines->InsertNextCell(line0);
	}

	// Add the lines to the polydata container
	lines_poly_data->SetLines(lines);
	line_mapper->Update();
	//hint_dialog->setHint(-1);
	refreshVisibility();
	emit cancelHint();
	emit cancelHint3D();
}

void VTKViewer::getCenterPointFromCellId(vtkIdType cell_id, double* p)
{
	vtkPoints* pts1 = poly_data->GetCell(cell_id)->GetPoints();
	p[0] = p[1] = p[2] = 0.0f;
	for (int i = 0; i < 3; i++)
	{
		double cords1[3];
		pts1->GetPoint(i, cords1);
		p[0] += cords1[0];
		p[1] += cords1[1];
		p[2] += cords1[2];
	}
	p[0] /= 3;
	p[1] /= 3;
	p[2] /= 3;
}

/*
void VTKViewer::enableCancelButton()
{
	cancel_button->setEnabled(true);
}

void VTKViewer::disableCancelButton()
{
	cancel_button->setEnabled(false);
}
*/

std::vector<cv::Point3d> VTKViewer::updateMainWindow()
{
	// return a pts_conf_3d
	std::vector<cv::Point3d> result;
	for (int i = 0; i < TOTAL_POINTS_IN_SINGLE_VIEW; i++)
	{
		double p[3];
		getCenterPointFromCellId(all_landmarks_cell_id->GetValue(i), p);
		result.push_back(cv::Point3d(p[0], p[1], p[2])); // what is the default confidence?
	}
	return result;
}

/*
void VTKViewer::getAllVisiblePointsInOneView()
{
	// set a new output port



	// hide this output port

	// set camera matrix

	// get visible points

	// return
}
*/

/*
void VTKViewer::countPoints()
{
	this->qvtkWidget->update();
	this->qvtkWidget->GetRenderWindow()->Render();
	this->visible_filter->Update();
	qDebug() << "There are currently: " << this->visible_filter->GetOutput()->GetNumberOfPoints() << " visible.";
	qDebug() << this->visible_filter->GetOutput()->GetNumberOfCells();
}
*/

/*
void VTKViewer::countPointsUsingHardwareSelector()
{
	vtkSmartPointer<vtkHardwareSelector> selector = vtkSmartPointer<vtkHardwareSelector>::New();
	selector->SetRenderer(renderer);
	int* temp = qvtkWidget->GetRenderWindow()->GetSize();
	unsigned int windowSize[4];
	windowSize[0] = temp[2];
	windowSize[1] = temp[3];
	windowSize[2] = temp[0];
	windowSize[3] = temp[1];
	
	for(unsigned int i = 0; i < 4; i++)
	{
	windowSize[i] = temp[i];
	}
	
	selector->SetArea(windowSize);
	selector->SetFieldAssociation(vtkDataObject::FIELD_ASSOCIATION_CELLS);
	vtkSelection* selection = selector->Select();
	qDebug() << "Selection has " << selection->GetNumberOfNodes() << " nodes.";

	vtkSmartPointer<vtkExtractSelection> extractSelection =
		vtkSmartPointer<vtkExtractSelection>::New();

	extractSelection->SetInputData(0, poly_data);
	extractSelection->SetInputData(1, selection);
	extractSelection->Update();
	// GetNumberOfElements: 1: cells
	// GetNumberOfElements: 0: points
	qDebug() << extract_selection->GetOutputDataObject(0)->GetNumberOfElements(1);



	
	vtkSmartPointer<vtkDataSetMapper> mapper =
		vtkSmartPointer<vtkDataSetMapper>::New();
	mapper->SetInputConnection(extractSelection->GetOutputPort());

	
	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(1, 0, 0);
	this->Renderer->AddActor(actor);
	
}
*/

void VTKViewer::enableFineAnnotationButton()
{
	correct_this_point_button->setEnabled(true);
}

void VTKViewer::correctThisPointButtonPressed_()
{
	emit correctThisPointButtonPressed(now_chosen_id_in_array);
}

void VTKViewer::to2DButtonPressed_()
{
	emit to2DButtonPressed();
}