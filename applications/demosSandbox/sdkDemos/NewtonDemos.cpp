/* Copyright (c) <2009> <Newton Game Dynamics>
* 
* This software is provided 'as-is', without any express or implied
* warranty. In no event will the authors be held liable for any damages
* arising from the use of this software.
* 
* Permission is granted to anyone to use this software for any purpose,
* including commercial applications, and to alter it and redistribute it
* freely
*/

// NewtonDemos.cpp : Defines the entry point for the application.
//


#include <toolbox_stdafx.h>
#include "SkyBox.h"
#include "DemoCamera.h"
#include "NewtonDemos.h"
#include "PhysicsUtils.h"
#include "DebugDisplay.h"
#include "DemoEntityManager.h"



#define DEFAULT_SCENE	0			// using NetwonMesh Tool
//#define DEFAULT_SCENE	1			// Coefficients of friction
//#define DEFAULT_SCENE	2			// Coefficients of restitution
//#define DEFAULT_SCENE	3			// Precessing tops
//#define DEFAULT_SCENE	4			// closest distance
//#define DEFAULT_SCENE	5			// primitive collision
//#define DEFAULT_SCENE	6 			// Kinematic bodies
//#define DEFAULT_SCENE	7			// primitive convex cast 
//#define DEFAULT_SCENE	8			// Box stacks
//#define DEFAULT_SCENE	9			// simple level mesh collision
//#define DEFAULT_SCENE	10			// optimized level mesh collision
//#define DEFAULT_SCENE	11			// height field Collision
//#define DEFAULT_SCENE	12			// infinite user plane collision
//#define DEFAULT_SCENE	13			// user height field Collision
//#define DEFAULT_SCENE	14			// compound Collision
//#define DEFAULT_SCENE	15			// pjani compound bug
//#define DEFAULT_SCENE	16			// uniform Scaled Collision
//#define DEFAULT_SCENE	17			// non Uniform Scaled Collision
//#define DEFAULT_SCENE	18			// scaled mesh collision
//#define DEFAULT_SCENE	19			// simple convex decomposition
//#define DEFAULT_SCENE	20			// scene Collision
//#define DEFAULT_SCENE	21          // simple boolean operators 
//#define DEFAULT_SCENE	22			// simple convex Shatter
//#define DEFAULT_SCENE	23			// multi ray casting using the threading Job scheduler
//#define DEFAULT_SCENE	24			// continue collision
//#define DEFAULT_SCENE	25			// puck slide continue collision
//#define DEFAULT_SCENE	26			// basic rag doll
//#define DEFAULT_SCENE	27			// basic car
//#define DEFAULT_SCENE	28			// high performance super car
//#define DEFAULT_SCENE	29			// basic player controller
//#define DEFAULT_SCENE	30			// advanced player controller
//#define DEFAULT_SCENE	31			// cloth patch			
//#define DEFAULT_SCENE	32			// soft bodies			


void Friction (DemoEntityManager* const scene);
void Restitution (DemoEntityManager* const scene);
void PrecessingTops (DemoEntityManager* const scene);
void ClosestDistance (DemoEntityManager* const scene);
void ConvexCast (DemoEntityManager* const scene);
void PrimitiveCollision (DemoEntityManager* const scene);
void KinematicBodies (DemoEntityManager* const scene);
void ClothPath(DemoEntityManager* const scene);
void SoftBodies (DemoEntityManager* const scene);
void BasicBoxStacks (DemoEntityManager* const scene);
void SimpleMeshLevelCollision (DemoEntityManager* const scene);
void OptimizedMeshLevelCollision (DemoEntityManager* const scene);
void UniformScaledCollision (DemoEntityManager* const scene);
void NonUniformScaledCollision (DemoEntityManager* const scene);
void ScaledMeshCollision (DemoEntityManager* const scene);
void ContinueCollision (DemoEntityManager* const scene);
void PuckSlide (DemoEntityManager* const scene);
void SceneCollision (DemoEntityManager* const scene);
void CompoundCollision(DemoEntityManager* const scene);
void PostCompoundCreateBuildTest(DemoEntityManager* const scene);
void SimpleConvexApproximation(DemoEntityManager* const scene);
void SimpleBooleanOperations(DemoEntityManager* const scene);
void SimpleConvexShatter (DemoEntityManager* const scene);
void UsingNewtonMeshTool (DemoEntityManager* const scene);
void MultiRayCast (DemoEntityManager* const scene);
void BasicCar (DemoEntityManager* const scene);
void SuperCar (DemoEntityManager* const scene);
void BasicPlayerController (DemoEntityManager* const scene);
void AdvancedPlayerController (DemoEntityManager* const scene);
void HeightFieldCollision (DemoEntityManager* const scene);
void UserPlaneCollision (DemoEntityManager* const scene);
void UserHeightFieldCollision (DemoEntityManager* const scene);
void DescreteRagDoll (DemoEntityManager* const scene);


NewtonDemos::SDKDemos NewtonDemos::m_demosSelection[] = 
{
	{wxT("Using the newton mesh tool"), wxT("demonstrate how to use the newton mesh toll for mesh manipulation"), UsingNewtonMeshTool},
	{wxT("Coefficients of friction"), wxT("demonstrate the effect of various coefficient of friction"), Friction},
	{wxT("Coefficients of restitution"), wxT("demonstrate the effect of various coefficient of restitution"), Restitution},
	{wxT("Precessing tops"), wxT("show natural precession"), PrecessingTops},
	{wxT("Closest distance"), wxT("demonstrate closest distance to a convex shape"), ClosestDistance},
	{wxT("Primitive Collision"), wxT("demonstrate separate collision of primitives"), PrimitiveCollision},
	{wxT("Kinematic bodies"), wxT("demonstrate separate collision of primitives"), KinematicBodies},
	{wxT("Primitive convex cast"), wxT("demonstrate separate primitive convex cast"), ConvexCast},
	{wxT("Simple box Stacks"), wxT("show simple stack of Boxes"), BasicBoxStacks},
	{wxT("Unoptimized mesh collision"), wxT("show simple level mesh"), SimpleMeshLevelCollision},
	{wxT("Optimized mesh collision"), wxT("show optimized level mesh"), OptimizedMeshLevelCollision},
	{wxT("Height field collision mesh"), wxT("show high file collision mesh"), HeightFieldCollision},
	{wxT("User infinite Plane collision mesh"), wxT("show high file collision mesh"), UserPlaneCollision},
	{wxT("User Height field collision mesh"), wxT("show high file collision mesh"), UserHeightFieldCollision},
	{wxT("Compound collision shape"), wxT("demonstrate compound collision"), CompoundCollision},
	{wxT("PostCompoundCreateBuildTest"), wxT("PostCompoundCreateBuildTest"), PostCompoundCreateBuildTest},
	{wxT("Uniform scaled collision shape"), wxT("demonstrate scaling shape"), UniformScaledCollision},
	{wxT("Non uniform scaled collision shape"), wxT("demonstrate scaling shape"), NonUniformScaledCollision},
	{wxT("Scaled mesh collision"), wxT("demonstrate scaling mesh scaling collision"), ScaledMeshCollision},
	{wxT("Simple convex decomposition"), wxT("demonstrate convex decomposition and compound collision"), SimpleConvexApproximation},
	{wxT("Multi geometry collision"), wxT("show static mesh with the ability of moving internal parts"), SceneCollision},
	{wxT("Simple boolean operations"), wxT("demonstrate simple boolean operations "), SimpleBooleanOperations},
	{wxT("Simple convex Shatter"), wxT("demonstrate fracture destruction using Voronoi partition"), SimpleConvexShatter},
	{wxT("Parallel ray cast"), wxT("using the threading Job scheduler"), MultiRayCast},
	{wxT("Continue collision"), wxT("show continue collision"), ContinueCollision},
	{wxT("Puck slide"), wxT("show continue collision"), PuckSlide},
	{wxT("Basic ragdoll"), wxT("demonstrate simple rag doll"), DescreteRagDoll},
	{wxT("Basic car"), wxT("implement a basic car"), BasicCar},
	//{wxT("High performance super car"), wxT("implement a high performance ray cast car"), SuperCar},
	{wxT("High performance super car"), wxT("implement a high performance ray cast car"), BasicCar},
	{wxT("Basic player controller"), wxT("demonstrate simple player controller"), BasicPlayerController},
	{wxT("Advanced player controller"), wxT("demonstrate player interacting with other objects"), AdvancedPlayerController},
	{wxT("Simple cloth Path"), wxT("show simple cloth path"), ClothPath},
	{wxT("Simple soft Body"), wxT("show simple soft body"), SoftBodies},

//	{wxT("basic convex hull stacking"), wxT("demonstrate convex hull stacking"), BasicConvexStacks},
//	{wxT("basic unstable stacking"), wxT("demonstrate stability stacking unstable objects"), UnstableStacks},
//	{wxT("Jenga stacking"), wxT("demonstrate Jenga game"), Jenga},
//	{wxT("Large Jenga stacking"), wxT("demonstrate Jenga game"), JengaTall},
//	{wxT("small pyramid stacking"), wxT("demonstrate small pyramid stacking"), CreatePyramid},
//	{wxT("wall stacking"), wxT("demonstrate wall stacking"), CreateWalls},
//	{wxT("small tower stacking"), wxT("demonstrate tower stacking"), CreateTower},
//	{wxT("large tower stacking"), wxT("demonstrate tower stacking"), CreateTowerTall},
//	{wxT("user defined polygon static collision"), wxT("demonstrate user defined polygon static collision"), UserHeighMapColliion},
//	{wxT("attractive magnets force field"), wxT("demonstrate attractive force field"), Magnets},
//	{wxT("repulsive magnets force field"), wxT("demonstrate repulsive magnet force field"), Repulsive},
//	{wxT("Archimedes buoyancy force field"), wxT("demonstrate user define Archimedes as force field"), ArchimedesBuoyancy},
//	{wxT("legacy joints"), wxT("demonstrate the build in joints"), LegacyJoints},
//	{wxT("custom joints"), wxT("demonstrate custom joints"), BasicCustomJoints},
//	{wxT("Simple robots"), wxT("demonstrate custom joints robot"), BasicRobots},
//	{wxT("motorized robots"), wxT("demonstrate motorized custom joints robot"), TracktionJoints},

//	{wxT("skinned rag doll"), wxT("demonstrate simple rag doll"), SkinRagDoll},
};


int NewtonDemos::m_threadsTracks[] = {1, 2, 3, 4, 8, 12, 16};



class NewtonDemosApp: public wxApp
{
	virtual bool OnInit()
	{
		// check for memory leaks
		#if defined(_DEBUG) && defined(_MSC_VER)
			// Track all memory leaks at the operating system level.
			// make sure no Newton tool or utility leaves leaks behind.
			_CrtSetDbgFlag(_CRTDBG_LEAK_CHECK_DF|_CrtSetDbgFlag(_CRTDBG_LEAK_CHECK_DF));
		#endif

		// Set the memory allocation function before creation the newton world
		// this is the only function that can be called before the creation of the newton world.
		// it should be called once, and the the call is optional 
		NewtonSetMemorySystem (PhysicsAlloc, PhysicsFree);

		int version = NewtonWorldGetVersion();
		wxString tittle;
		tittle.Printf (wxT ("Newton %d.%02d SDK demos"), version / 100, version % 100);
		NewtonDemos* const frame = new NewtonDemos(tittle, wxDefaultPosition, wxSize(1024, 768));
		

		frame->Show(true);
		SetTopWindow(frame);

		// initialize open gl graphics
		if (frame->m_scene) {
			frame->m_scene->InitGraphicsSystem();
		}

		// load the default Scene		
		//frame->LoadDemo (DEFAULT_SCENE);
		wxMenuEvent loadDemo (wxEVT_COMMAND_MENU_SELECTED, NewtonDemos::ID_RUN_DEMO + DEFAULT_SCENE);
		frame->GetEventHandler()->ProcessEvent(loadDemo);

		return true;
	}

	// memory allocation for Newton
	static void* PhysicsAlloc (int sizeInBytes)
	{
		m_totalMemoryUsed += sizeInBytes;
		return new char[sizeInBytes];
	}

	// memory free use by the engine
	static void PhysicsFree (void* ptr, int sizeInBytes)
	{
		m_totalMemoryUsed -= sizeInBytes;
		delete[] (char*)ptr;
	}

	static int m_totalMemoryUsed;
};

int NewtonDemosApp::m_totalMemoryUsed = 0;


IMPLEMENT_APP(NewtonDemosApp)

BEGIN_EVENT_TABLE(NewtonDemos, wxFrame)
	// mandatory menu events for mac cocoa osx  support

	EVT_MENU(wxID_ABOUT, NewtonDemos::OnAbout)
	EVT_MENU(wxID_EXIT, NewtonDemos::OnQuit)
	EVT_MENU(wxID_HELP, NewtonDemos::OnAbout)
	EVT_MENU(wxID_PREFERENCES, NewtonDemos::OnAbout)
	EVT_MENU(wxID_NEW, NewtonDemos::OnNew)

	// game menus events
	EVT_MENU_RANGE(ID_RUN_DEMO, ID_RUN_DEMO_RANGE, NewtonDemos::OnRunDemo)

	EVT_MENU(ID_AUTOSLEEP_MODE,	NewtonDemos::OnAutoSleepMode)
	EVT_MENU(ID_SHOW_STATISTICS, NewtonDemos::OnShowStatistics)
	EVT_MENU(ID_USE_PARALLEL_SOLVER, NewtonDemos::OnUseParallelSolver)

	EVT_MENU(ID_HIDE_VISUAL_MESHES,	NewtonDemos::OnHideVisualMeshes)

	EVT_MENU_RANGE(ID_SHOW_COLLISION_MESH, ID_SHOW_COLLISION_MESH_RANGE, NewtonDemos::OnShowCollisionLines)

	EVT_MENU(ID_SHOW_CONTACT_POINTS, NewtonDemos::OnShowContactPoints)
	EVT_MENU(ID_SHOW_NORMAL_FORCES,	NewtonDemos::OnShowNormalForces)
	EVT_MENU(ID_SHOW_AABB, NewtonDemos::OnShowAABB)
	EVT_MENU(ID_SHOW_CENTER_OF_MASS, NewtonDemos::OnShowCenterOfMass)
	EVT_MENU(ID_SHOW_JOINTS, NewtonDemos::OnShowShowJoints)
	

	EVT_MENU_RANGE(ID_PLATFORMS, ID_PLATFORMS_MAX, NewtonDemos::OnSelectHardwareDevice)

	EVT_MENU(ID_SHOW_CONCURRENCE_PROFILER, NewtonDemos::OnShowConcurrentProfiler)
	EVT_MENU(ID_SHOW_PROFILER,	NewtonDemos::OnShowThreadProfiler)

	EVT_MENU(ID_SELECT_ALL_PROFILERS, NewtonDemos::OnSelectAllPerformanceChart)
	EVT_MENU(ID_UNSELECT_ALL_PROFILERS,	NewtonDemos::OnUnselectAllPerformanceChart)
	EVT_MENU_RANGE (ID_SHOW_PHYSICS_PROFILER, ID_SHOW_PHYSICS_PROFILER_COUNT, NewtonDemos::OnShowProfiler)

	EVT_MENU(ID_CONCURRENT_PHYSICS_UPDATE, NewtonDemos::OnRunPhysicsConcurrent)
	EVT_MENU_RANGE(ID_SELECT_MICROTHREADS, ID_SELECT_MICROTHREADS_COUNT, NewtonDemos::OnSelectNumberOfMicroThreads)

	
	EVT_MENU(ID_SERIALIZE, NewtonDemos::OnSerializeWorld)
	EVT_MENU(ID_DESERIALIZE, NewtonDemos::OnDeserializeWorld)



//	FXMAPFUNC(SEL_COMMAND,		NewtonDemos::ID_LOAD,								NewtonDemos::onLoad),
//	FXMAPFUNC(SEL_COMMAND,		NewtonDemos::ID_SAVE,								NewtonDemos::onSave),


END_EVENT_TABLE()


NewtonDemos::NewtonDemos(const wxString& title, const wxPoint& pos, const wxSize& size)
	:wxFrame(NULL, -1, title, pos, size)
	,m_mainMenu(NULL)
	,m_statusbar(NULL)
	,m_scene(NULL)
	,m_physicsUpdateMode(0)
	,m_suspendVisualUpdates(true)
	,m_autoSleepState(true)
	,m_useParallelSolver(false)
	,m_hideVisualMeshes(false)
	,m_showContactPoints(false)
	,m_showNormalForces(false)
	,m_showAABB(false)
	,m_showJoints(false)
	,m_showCenterOfMass(false)
	,m_showStatistics(false)
	,m_concurrentProfilerState(false)
	,m_threadProfilerState(false)
	,m_hasJoysticController(false)
	,m_debugDisplayMode(0)
	,m_mousePosX(0)
	,m_mousePosY(0)
	,m_joytickX(0)
	,m_joytickY(0)
	,m_joytickButtonMask(0)
	,m_framesCount(0)
	,m_microthreadIndex(0)
	,m_hardwareDevice(0)
	,m_timestepAcc(0)
	,m_fps(0.0f)
{
	memset (m_profilerTracksMenu, 0, sizeof (m_profilerTracksMenu));

	// clear the key map
	memset (m_key, 0, sizeof (m_key));
	for (int i = 0; i < int (sizeof (m_keyMap)/sizeof (m_keyMap[0])); i ++) {
		m_keyMap[i] = i;
	}
	for (int i = 'a'; i <= 'z'; i ++) {
		m_keyMap[i] = i - 'a' + 'A';
	}

#ifdef WIN32
	m_keyMap[0] = VK_LBUTTON;
	m_keyMap[1] = VK_RBUTTON;
	m_keyMap[2] = VK_MBUTTON; 
#endif

	m_scene = new DemoEntityManager(this);
	m_statusbar = CreateStatusBar();
	int widths[] = {150, 160, 150, 90, 80, 100, 100};
	m_statusbar->SetFieldsCount (sizeof (widths)/sizeof (widths[0]), widths);
	CalculateFPS(0.0f);
	m_mainMenu = CreateMainMenu();

m_showStatistics = true;
//m_debugDisplayMode = 2;
//m_autoSleepState = false;
//m_useParallelSolver = true;
//m_scene->m_showProfiler[6] = 1;
m_scene->m_showProfiler[0] = 1;
}


NewtonDemos::~NewtonDemos()
{
}


wxMenuBar* NewtonDemos::CreateMainMenu()
{
	wxMenuBar* const mainMenu =  new wxMenuBar();

	// adding the file menu
	{
		wxMenu* const fileMenu = new wxMenu;

		fileMenu->Append(wxID_ABOUT, wxT("About"));
		
		fileMenu->AppendSeparator();
		fileMenu->Append(wxID_PREFERENCES, wxT("Preferences"));

		fileMenu->AppendSeparator();
		fileMenu->Append(wxID_NEW, wxT("&New"), wxT("Create a blank new scene"));

		fileMenu->AppendSeparator();
		fileMenu->Append(wxID_OPEN, wxT("&Open"), wxT("Open visual scene in dScene newton format"));
		fileMenu->Append(wxID_SAVE, wxT("&Save"), wxT("Save visual scene in dScene newton format"));

		fileMenu->AppendSeparator();
		fileMenu->Append(ID_SERIALIZE, wxT("&Serialize"), wxT("Serialize scene to binary file"));
		fileMenu->Append(ID_DESERIALIZE, wxT("&Deserialize"), wxT("Load previuoslly serialized scame"));

	//	fileMenu->AppendSeparator();
	//	fileMenu->Append(m_idImportPhysics, wxT("&Open physics scene"), wxT("Open physics scene in collada format"));
	//	fileMenu->Append(m_idExportPhysics, wxT("&Save physics scene"), wxT("Save physics in collada format"));

		fileMenu->AppendSeparator();
		fileMenu->Append(wxID_EXIT, wxT("E&xit\tAlt-X"), wxT("Quit SDK sample") );

		// add main menus to menu bar
		mainMenu->Append(fileMenu, wxT("&File"));
	}

	// engine all demo examples
	{
		wxMenu* const sdkDemos = new wxMenu;
		int demosCount = int (sizeof (m_demosSelection) / sizeof m_demosSelection[0]);
		for (int i = 0; i < demosCount; i ++) {
			sdkDemos->AppendRadioItem (NewtonDemos::ID_RUN_DEMO + i,  m_demosSelection[i].m_name, m_demosSelection[i].m_description);
		}

		mainMenu->Append(sdkDemos, wxT("&Demos"));
	}

	// option menu
	{
		wxMenu* const optionsMenu = new wxMenu;;

		optionsMenu->AppendCheckItem(ID_AUTOSLEEP_MODE, wxT("Auto sleep mode"), wxT("toogle auto sleep bodies"));
		optionsMenu->Check (ID_AUTOSLEEP_MODE, m_autoSleepState);

		optionsMenu->AppendCheckItem(ID_SHOW_STATISTICS, wxT("Show Stats on screen"), wxT("toogle on screen frame rate and other stats"));
		optionsMenu->AppendCheckItem(ID_USE_PARALLEL_SOLVER, wxT("Parallel solver on"));

		optionsMenu->AppendSeparator();
		optionsMenu->AppendRadioItem(ID_SHOW_COLLISION_MESH, wxT("Hide collision Mesh"));
		optionsMenu->AppendRadioItem(ID_SHOW_COLLISION_MESH + 1, wxT("Show solid collision Mesh"));
		optionsMenu->AppendRadioItem(ID_SHOW_COLLISION_MESH + 2, wxT("Show wire frame collision Mesh"));

		optionsMenu->AppendSeparator();
		optionsMenu->AppendCheckItem(ID_HIDE_VISUAL_MESHES, wxT("Hide visual meshes"));
		optionsMenu->AppendCheckItem(ID_SHOW_CONTACT_POINTS, wxT("Show contact points"));
		optionsMenu->AppendCheckItem(ID_SHOW_NORMAL_FORCES, wxT("Show normal forces"));
		optionsMenu->AppendCheckItem(ID_SHOW_AABB, wxT("Show aabb"));
		optionsMenu->AppendCheckItem(ID_SHOW_CENTER_OF_MASS, wxT("Show center of mass"));
		optionsMenu->AppendCheckItem(ID_SHOW_JOINTS, wxT("show Joint debug info"));
	

		optionsMenu->AppendSeparator();
		int platformsCount = NewtonEnumrateDevices (m_scene->GetNewton());
		for (int i = 0; i < platformsCount; i ++) {
			wxString label;
			char platform[256];
			

			NewtonGetDeviceString (m_scene->GetNewton(), i, platform, sizeof (platform));
			#ifdef _POSIX_VER
				wxChar wPlatform[256];
				mbstowcs (wPlatform, platform, sizeof (platform));
				wxString tmp (wPlatform);
				label.Printf (wxT(" hardware mode %s"), tmp.c_str());
			#else 
				label.Printf (wxT(" hardware mode %s"), wxString(platform));
			#endif
			optionsMenu->AppendRadioItem(ID_PLATFORMS + i, label);
		}
		//optionsMenu->Check(ID_PLATFORMS, true);

		optionsMenu->AppendSeparator();
		optionsMenu->Append(ID_SHOW_CONCURRENCE_PROFILER, wxT("Show concurrent profiler"));
		optionsMenu->Append(ID_SHOW_PROFILER, wxT("Show micro thread profiler"));

		optionsMenu->AppendSeparator();
		optionsMenu->Append(ID_SELECT_ALL_PROFILERS, wxT("select all profiler"));
		optionsMenu->Append(ID_UNSELECT_ALL_PROFILERS, wxT("unselect all profiler"));

		wxMenu* const profilerSubMenu = new wxMenu;
		m_profilerTracksMenu[0] = profilerSubMenu->AppendCheckItem(ID_SHOW_PHYSICS_PROFILER + 0, wxT("show global physics update performance chart"));
		m_profilerTracksMenu[1] = profilerSubMenu->AppendCheckItem(ID_SHOW_PHYSICS_PROFILER + 1, wxT("global collision update performance chart"));
		m_profilerTracksMenu[2] = profilerSubMenu->AppendCheckItem(ID_SHOW_PHYSICS_PROFILER + 2, wxT("broad phase collision performance chart"));
		m_profilerTracksMenu[3] = profilerSubMenu->AppendCheckItem(ID_SHOW_PHYSICS_PROFILER + 3, wxT("narrow phase collision performance chart"));
		m_profilerTracksMenu[4] = profilerSubMenu->AppendCheckItem(ID_SHOW_PHYSICS_PROFILER + 4, wxT("global dynamics update performance chart"));
		m_profilerTracksMenu[5] = profilerSubMenu->AppendCheckItem(ID_SHOW_PHYSICS_PROFILER + 5, wxT("dynamics setup performance chart"));
		m_profilerTracksMenu[6] = profilerSubMenu->AppendCheckItem(ID_SHOW_PHYSICS_PROFILER + 6, wxT("dynamics solver performance chart"));
		m_profilerTracksMenu[7] = profilerSubMenu->AppendCheckItem(ID_SHOW_PHYSICS_PROFILER + 7, wxT("force and torque callback performance chart"));
		m_profilerTracksMenu[8] = profilerSubMenu->AppendCheckItem(ID_SHOW_PHYSICS_PROFILER + 8, wxT("pre-simulation listener"));
		m_profilerTracksMenu[9] = profilerSubMenu->AppendCheckItem(ID_SHOW_PHYSICS_PROFILER + 9, wxT("post-simulation listener"));
//		if (mainFrame->m_physicProfilerState) {
//			m_profilerTracksMenu[0]->setCheck(true);
//		}
		optionsMenu->AppendSubMenu (profilerSubMenu, wxT("select sub profiler"));


		optionsMenu->AppendSeparator();
		optionsMenu->AppendCheckItem(ID_CONCURRENT_PHYSICS_UPDATE, wxT("Concurrent physics update"));

		wxMenu* const microThreadedsSubMenu = new wxMenu;
		for (int i = 0 ; i < int (sizeof (m_threadsTracks)/ sizeof (m_threadsTracks[0])); i ++) {
			wxString msg;
			msg.Printf(wxT ("%d micro threads"), m_threadsTracks[i]);
			microThreadedsSubMenu->AppendRadioItem(ID_SELECT_MICROTHREADS + i, msg);
		}
		optionsMenu->AppendSubMenu (microThreadedsSubMenu, wxT("select microThread count"));


		mainMenu->Append(optionsMenu, wxT("&Options"));
	}

	// add help menu
	{
		wxMenu* const helpMenu = new wxMenu;;

		helpMenu->Append(wxID_HELP, wxT("About"));
//		helpMenu->Append(NewtonDemos::ID_ON_ABOUT, wxT("About"));
		mainMenu->Append(helpMenu, wxT("&Help"));
	}

	SetMenuBar(mainMenu);
	return mainMenu;
}


void NewtonDemos::BEGIN_MENU_OPTION()																				
{
	m_suspendVisualUpdates = true;																			
	if (m_scene && m_scene->GetNewton()) {																			
		NewtonWaitForUpdateToFinish (m_scene->GetNewton());												
	}
}


void NewtonDemos::END_MENU_OPTION()
{
	m_suspendVisualUpdates = false;																			
	if (m_scene && m_scene->GetNewton()) {		
		NewtonWaitForUpdateToFinish (m_scene->GetNewton());
		SetAutoSleepMode (m_scene->GetNewton(), !m_autoSleepState);
		NewtonSetMultiThreadSolverOnSingleIsland (m_scene->GetNewton(), m_useParallelSolver ? 1 : 0);	
	}
}


void NewtonDemos::RestoreSettings ()
{
	NewtonSetCurrentDevice (m_scene->GetNewton(), m_hardwareDevice); 
	NewtonSetThreadsCount(m_scene->GetNewton(), m_threadsTracks[m_microthreadIndex]); 
}


void NewtonDemos::LoadDemo (int index)
{
	BEGIN_MENU_OPTION();

	m_scene->Cleanup();

	m_scene->SetCurrent();
	m_demosSelection[index].m_launchDemoCallback (m_scene);
	m_scene->SwapBuffers(); 

	RestoreSettings ();
	m_scene->ResetTimer();

	END_MENU_OPTION();
}



void NewtonDemos::KeyDown(const wxKeyEvent &event)
{
	int keyCode = event.GetKeyCode();
	if (keyCode == WXK_ESCAPE)  {
		// send a display refresh event in case the runtime update is stopped bu the user.
		wxMenuEvent exitEvent (wxEVT_COMMAND_MENU_SELECTED, wxID_EXIT);
		GetEventHandler()->ProcessEvent(exitEvent);
	}


	if (!event.GetModifiers()) {
		int code = keyCode & 0xff; 
		m_key[m_keyMap[code]] = true;
	}
}


void NewtonDemos::KeyUp(const wxKeyEvent &event)
{
	if (!event.GetModifiers()) {
		int keyCode = event.GetKeyCode();
		int code = keyCode & 0xff;
		m_key[m_keyMap[code]] = false;
	}
}


void NewtonDemos::MouseAction(const wxMouseEvent &event)
{
	m_mousePosX = event.GetX();
	m_mousePosY = event.GetY();

	if (event.LeftIsDown()) {
		m_key[m_keyMap[0]] = true;
	} else {
		m_key[m_keyMap[0]] = false;
	}

	if (event.RightIsDown()) {
		m_key[m_keyMap[1]] = true;
	} else {
		m_key[m_keyMap[1]] = false;
	}
}

bool NewtonDemos::GetKeyState( int key) const
{
	return m_key[m_keyMap[key & 0xff]] ? true : false;
}

bool NewtonDemos::GetMousePosition (int& posX, int& posY) const
{
	posX = m_mousePosX;
	posY = m_mousePosY;
	return true;
}

bool NewtonDemos::GetJoytickPosition (dFloat& posX, dFloat& posY, int& buttonsMask) const
{
	buttonsMask = m_joytickButtonMask;
	posX = dFloat (m_joytickX - 32767) / 32768.0f;
	posY = -dFloat (m_joytickY - 32767) / 32768.0f;
	return m_hasJoysticController;
}


bool NewtonDemos::GetMouseKeyState (int button) const
{
	if ((button >= 0) && (button <= 2)) {
		return m_key[m_keyMap[button]] ? true : false;
	}

	return false;
}

void NewtonDemos::CalculateFPS(float timestep)
{
	m_framesCount ++;
	m_timestepAcc += timestep;

	// this probably happing on loading of and a pause, just rest counters
	if ((m_timestepAcc <= 0.0f) || (m_timestepAcc > 2.0f)){
		m_timestepAcc = 0;
		m_framesCount = 0;
	}

	//update fps every quarter of a second
	if (m_timestepAcc >= 0.25f) {
		m_fps = float (m_framesCount) / m_timestepAcc;
		m_timestepAcc -= 0.25f;
		m_framesCount = 0.0f;

		wxString statusText;
		NewtonWorld* const world = m_scene->GetNewton();
		char platform[256];
		NewtonGetDeviceString(world, NewtonGetCurrentDevice(world), platform, sizeof (platform));
		int memoryUsed = NewtonGetMemoryUsed() / (1024) ;
		
		statusText.Printf (wxT ("render fps: %7.2f"), m_fps);
		m_statusbar->SetStatusText (statusText, 0);

		statusText.Printf (wxT ("physics time: %4.2f ms"), m_scene->GetPhysicsTime() * 1000.0f);
		m_statusbar->SetStatusText (statusText, 1);

		statusText.Printf (wxT ("memory: %d kbytes"), memoryUsed);
		m_statusbar->SetStatusText (statusText, 2);

		statusText.Printf (wxT ("bodies: %d"), NewtonWorldGetBodyCount(world));
		m_statusbar->SetStatusText (statusText, 3);

		statusText.Printf (wxT ("threads: %d"), NewtonGetThreadsCount(world));
		m_statusbar->SetStatusText (statusText, 4);

		statusText.Printf (wxT ("auto sleep: %s"), m_autoSleepState ? "on" : "off");
		m_statusbar->SetStatusText (statusText, 5);

		char floatMode[128];
		NewtonGetDeviceString (m_scene->GetNewton(), m_hardwareDevice, floatMode, sizeof (floatMode));
		statusText.Printf (wxT ("instructions: %s"), floatMode);
		m_statusbar->SetStatusText (statusText, 6);

	}
}




void NewtonDemos::OnAbout(wxCommandEvent& event)
{
	wxString msg;
	int version = NewtonWorldGetVersion();
	msg.Printf(wxT ("Hello to Newton Dynamics SDK %d.%02d"), version / 100, version % 100);
	wxMessageBox(msg, wxT ("Newton Dynanics"), wxOK | wxICON_INFORMATION, this);
}


void NewtonDemos::OnQuit(wxCommandEvent& event)
{
	Close ();
}

void NewtonDemos::OnRunDemo(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();

	wxMenu* const menu = m_mainMenu->GetMenu(m_mainMenu->FindMenu(_ ("Demos")));
	if (!menu->IsChecked(event.GetId())) {
		menu->Check (event.GetId(), true);
	}

	LoadDemo (event.GetId() - ID_RUN_DEMO);

	RestoreSettings ();
	END_MENU_OPTION();
}


void NewtonDemos::OnAutoSleepMode(wxCommandEvent& event)	
{
	BEGIN_MENU_OPTION();
	m_autoSleepState = event.IsChecked();
	END_MENU_OPTION();
}


void NewtonDemos::OnHideVisualMeshes(wxCommandEvent& event)	
{
	BEGIN_MENU_OPTION();
	m_hideVisualMeshes = event.IsChecked();
	END_MENU_OPTION();
}


void NewtonDemos::OnShowCollisionLines(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();
	m_debugDisplayMode = event.GetId() - ID_SHOW_COLLISION_MESH;
	SetDebugDisplayMode (m_debugDisplayMode);
	END_MENU_OPTION();
}


void NewtonDemos::OnShowContactPoints(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();
	m_showContactPoints = event.IsChecked(); 
	END_MENU_OPTION();
}

void NewtonDemos::OnShowNormalForces(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();
	m_showNormalForces = event.IsChecked(); 
	END_MENU_OPTION();
}

void NewtonDemos::OnShowAABB(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();
	m_showAABB = event.IsChecked(); 
	END_MENU_OPTION();
}

void NewtonDemos::OnShowCenterOfMass(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();
	m_showCenterOfMass = event.IsChecked(); 
	END_MENU_OPTION();
}

void NewtonDemos::OnShowShowJoints(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();
	m_showJoints = event.IsChecked(); 
	END_MENU_OPTION();
}


void NewtonDemos::OnUseParallelSolver(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();
	m_useParallelSolver = event.IsChecked(); 
	NewtonSetMultiThreadSolverOnSingleIsland (m_scene->GetNewton(), m_useParallelSolver ? 1 : 0);
	END_MENU_OPTION();
}


void NewtonDemos::OnShowConcurrentProfiler(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();
	m_concurrentProfilerState = event.IsChecked(); 
	END_MENU_OPTION();
}


void NewtonDemos::OnShowThreadProfiler(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();
	m_threadProfilerState = event.IsChecked(); 
	END_MENU_OPTION();
}


void NewtonDemos::OnShowProfiler(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();
	int track = event.GetId() - ID_SHOW_PHYSICS_PROFILER;
	int state = event.IsChecked(); 
	m_scene->m_showProfiler[track] = state;
	END_MENU_OPTION();
}


void NewtonDemos::OnSelectAllPerformanceChart(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();
	int count = sizeof (m_profilerTracksMenu) / sizeof (m_profilerTracksMenu[0]);
	for (int i = 0; i < count; i ++) {
		if (m_profilerTracksMenu[i]) {
			m_scene->m_showProfiler[i] = 1;
			m_profilerTracksMenu[i]->Check(true);
		}
	}
	END_MENU_OPTION();
}

void NewtonDemos::OnUnselectAllPerformanceChart(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();
	int count = sizeof (m_profilerTracksMenu) / sizeof (m_profilerTracksMenu[0]);
	for (int i = 0; i < count; i ++) {
		if (m_profilerTracksMenu[i]) {
			m_scene->m_showProfiler[i] = 0;
			m_profilerTracksMenu[i]->Check(false);
		}
	}
	END_MENU_OPTION();
}


void NewtonDemos::OnRunPhysicsConcurrent(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();
	m_physicsUpdateMode = event.IsChecked(); 
	END_MENU_OPTION();
}

void NewtonDemos::OnSelectNumberOfMicroThreads(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();
	m_microthreadIndex = event.GetId() - ID_SELECT_MICROTHREADS;
	NewtonSetThreadsCount(m_scene->GetNewton(), m_threadsTracks[m_microthreadIndex]);
	END_MENU_OPTION();
}

void NewtonDemos::OnSelectHardwareDevice(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();
	m_hardwareDevice = event.GetId() - ID_PLATFORMS;
	NewtonSetCurrentDevice (m_scene->GetNewton(), m_hardwareDevice);
	END_MENU_OPTION();
}

void NewtonDemos::OnShowStatistics(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();
	m_showStatistics = event.IsChecked(); 
	END_MENU_OPTION();
}

void NewtonDemos::OnNew (wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();
	m_scene->Cleanup();
	RestoreSettings ();
	m_scene->ResetTimer();
	END_MENU_OPTION();
}


void NewtonDemos::OnSerializeWorld (wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();

	wxFileDialog open (this, wxT("Export a Newton Dynamics Serialized Physics Scene"), wxT("../../../media"), wxT(""), wxT("*.bin"));
	if (open.ShowModal() == wxID_OK) {
		wxString currentDocPath (open.GetPath());
		m_scene->SerializedPhysicScene (currentDocPath.mb_str());
	}
	m_scene->ResetTimer();

	END_MENU_OPTION();
}

void NewtonDemos::OnDeserializeWorld(wxCommandEvent& event)
{
	BEGIN_MENU_OPTION();

	wxFileDialog save (this, wxT("Import a Newton Dynamics Serialized Physics Scene"), wxT("../../../media"), wxT(""), wxT("*.bin"));
	if (save.ShowModal() == wxID_OK) {
		wxString currentDocPath (save.GetPath());
		//m_scene->DeserializedPhysicScene (currentDocPath.c_str());
		m_scene->DeserializedPhysicScene (currentDocPath.mb_str());
		RestoreSettings ();
	}

	m_scene->ResetTimer();
	END_MENU_OPTION();
}


#if 0
long NewtonDemos::onLoad(FXObject* sender, FXSelector id, void* eventPtr)
{
	BEGIN_MENU_OPTION();

	const FXchar patterns[]="Newton Dynamics Files (*.ngd)";
	FXFileDialog open(this,"Load Newton Dynamics scene");
	open.setPatternList(patterns);
	open.setDirectory ("../../../media");
	if(open.execute()){

		m_scene->Cleanup();

		// load the scene from a ngd file format
		m_scene->makeCurrent();
		m_scene->LoadScene (open.getFilename().text());
		m_scene->makeNonCurrent();

		// add a sky box to the scene, make the first object
		m_scene->Addtop (new SkyBox());

		// place camera into position
		dMatrix camMatrix (GetIdentityMatrix());
		//		camMatrix.m_posit = dVector (-40.0f, 10.0f, 0.0f, 0.0f);
		camMatrix = dYawMatrix(-0.0f * 3.1416f / 180.0f);
		camMatrix.m_posit = dVector (-5.0f, 1.0f, -0.0f, 0.0f);
		m_scene->SetCameraMatrix(camMatrix, camMatrix.m_posit);

		RestoreSettings ();
	}


	m_scene->ResetTimer();
	END_MENU_OPTION();
	return 1;
}

long NewtonDemos::onSave(FXObject* sender, FXSelector id, void* eventPtr)
{
	return 1;
}

#endif