//------------------------------------------------------------------------------
// 1/2	:切換image/depth
// 3	:初始點
//------------------------------------------------------------------------------

#include "stdafx.h"
#include <strsafe.h>
#include <math.h>
#include <limits>
#include <Wincodec.h>
#include "resource.h"
#include <Kinect.h>
//-------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <math.h>
#include <memory>
#include <windows.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include "include/scene.h"
#include "DrawImg/imgmap.h"
#include "DrawImg/glsetting.h"
#include "lib3D/obj.h"
#include "FootTracking/tracking.h"
#include <thread>
#include <sstream>

//
//#pragma comment(linker, "/subsystem:\"windows\" /entry:\"mainCRTStartup\"")

//
#define _WINSOCK_DEPRECATED_NO_WARNINGS
#pragma comment(lib, "Ws2_32.lib")
#include <WinSock2.h>
//
#define window_width  1920
#define window_height 1080
// alignThreadNum:可整除1080,每個cpu core一個為佳,在此以八核心電腦扣掉兩條(mThread,mThread2),所以用5
#define alignThreadNum 1
// namespace
using namespace std;
using namespace FootTracker;
// Image
point3f *depth;
float zValue[1920 * 1080];
color3uc *img, *d_img, *imgTest;
unsigned int frame_count = 100;
unsigned int frame_num = 0;
// General variable
int draw = 2, mode = 7;
int testmode = 1;
int shoeNum = 0;
int adjustmode = 0;
Mesh *obj, *foot, *hide, *simple, *tempObj;
Mesh *objL, *footL, *hideL;
Mesh *obj1;
float m[16], glm[16], mL[16], glmL[16];
Tracker right_feet, left_feet;
Point3f_Array selectpoints, selectpointsL;
MeshBuffer mbuf, mbufL;
float icpe[100], icpeL[100];
thread alignThread[alignThreadNum];
bool needChangeShoe = false;
bool needChangeFoot = false;
string changeObjName = "";
string changeObjNameL = "";
string changeFootName = "";
string changeFootNameL = "";

int WinNumber = NULL;
float distancex = 0;				     //在平移矩陣(glTranslatef();)中使用 
float distancey = 0;
float distancez = 0;
float distancef = 0;
float calibX;
float calibY;
float calibZ;
float calibF;
float fovy;
float aspect;
int hintLX;
int hintLY;
int hintRX;
int hintRY;
Vector4 floor_Vector;
float viewval = 0;
float ratioval = 0;
float rate = 0;
float modelScale = 1.0;
int isinit = 0;
float minY = 0;
int FPS = 0;
bool initFlagR = false;
bool initFlagL = false;
int initTimer = 0;
// functions
void init();
void idle();
void main_loop();
void display();
void timer(int value);
void reshape(int w, int h);
void keyboard(unsigned char key, int x, int y);
void Mouse(int, int, int, int);

// occlusion improvement 
Point3f footanglepoint, footanglepoint2;

// data
IKinectSensor*		pSensor = nullptr;
IColorFrameReader*	pColorFrameReader = nullptr;
IDepthFrameReader*	pDepthFrameReader = nullptr;
ICoordinateMapper*	pCoordinateMapper = nullptr;
IBodyFrameReader*	pBodyFrameReader = nullptr;
int		iColorWidth = 0, iColorHeight = 0;
UINT	uDepthPointNum = 0;
UINT	uColorPointNum = 0;
UINT	uColorBufferSize = 0;

UINT16*	pDepthBuffer = nullptr;
BYTE*	pColorBuffer = nullptr;
CameraSpacePoint* pCSPoints = nullptr;
Vector4 floor_vector;

void AlignDepth(int threadIdx);
void AlignDepthMaps();
void ByteArraryToColorImg();
void RightFeetTracking();
void LeftFeetTracking();
void skt()
{
	char message[200];
	string tmp;
	int r;
	WSAData wsaData;
	WORD DLLVSERION;
	DLLVSERION = MAKEWORD(2, 1);//Winsocket-DLL 版本

	//用 WSAStartup 開始 Winsocket-DLL
	r = WSAStartup(DLLVSERION, &wsaData);

	//宣告 socket 位址資訊(不同的通訊,有不同的位址資訊,所以會有不同的資料結構存放這些位址資訊)
	SOCKADDR_IN addr;
	int addrlen = sizeof(addr);

	//建立 socket
	SOCKET sListen; //listening for an incoming connection
	SOCKET sConnect; //operating if a connection was found

	//AF_INET：表示建立的 socket 屬於 internet family
	//SOCK_STREAM：表示建立的 socket 是 connection-oriented socket 
	sConnect = socket(AF_INET, SOCK_STREAM, NULL);

	//設定位址資訊的資料
	addr.sin_addr.s_addr = inet_addr("127.0.0.1");
	addr.sin_family = AF_INET;
	addr.sin_port = htons(1234);

	//設定 Listen
	sListen = socket(AF_INET, SOCK_STREAM, NULL);
	::bind(sListen, (SOCKADDR*)&addr, sizeof(addr));
	listen(sListen, SOMAXCONN);//SOMAXCONN: listening without any limit

	//等待連線
	SOCKADDR_IN clinetAddr;
	while (true)
	{
		cout << "waiting..." << endl;
		sConnect = accept(sListen, (SOCKADDR*)&clinetAddr, &addrlen);

		if (sConnect != INVALID_SOCKET)
		{
			cout << "a connection was found" << endl;
			printf("server: got connection from %s\n", inet_ntoa(addr.sin_addr));

			ZeroMemory(message, 200);
			recv(sConnect, message, sizeof(message), 0);
			cout << message << endl;
			if (string(message).substr(0, 4) == "shoe")
			{
				int nameLength = string(message).length() - 4;
				string objName = string(message).substr(4, nameLength);
				std::cout << "Change model: " << objName.c_str() << std::endl;
				changeObjName = "data/Shoes/" + objName + "-R.obj";
				changeObjNameL = "data/Shoes/" + objName + "-L.obj";
				needChangeShoe = true;
			}
			else if (string(message).substr(0, 4) == "foot")
			{
				int nameLength = string(message).length() - 4;
				string objName = string(message).substr(4, nameLength);
				std::cout << "Change foot: " << objName.c_str() << std::endl;
				changeFootName = "data/" + objName + ".obj";
				changeFootNameL = "data/" + objName + "_left.obj";
				needChangeFoot = true;
			}
			else
			{
				cout << "ERROR" << endl;
			}

		}
	}
}
void SetWindowStyle()
{
	long dwStyle;
	HWND hwndGlut;
	hwndGlut = FindWindow(NULL, L"Virtual Try-on");
	// 設定無邊框視窗
	dwStyle = GetWindowLong(hwndGlut, GWL_STYLE);
	//dwStyle &= ~(WS_CAPTION | WS_THICKFRAME | WS_MINIMIZE | WS_MAXIMIZE | WS_SYSMENU);
	dwStyle &= ~(WS_CAPTION | WS_MINIMIZE | WS_MAXIMIZE | WS_SYSMENU);
	SetWindowLongPtr(hwndGlut, GWL_STYLE, dwStyle);

	SetWindowPos(hwndGlut, NULL, 0, 0, 0, 0, SWP_FRAMECHANGED |
		SWP_NOMOVE | SWP_NOSIZE | SWP_NOZORDER | SWP_NOOWNERZORDER);
}

void ReadINIFile()
{
	wchar_t value[MAX_PATH];
	value[0] = '\0';
	// [Model]
	GetPrivateProfileString(L"Model", L"ShoeName", L"data/yaya_shoe2.obj", value, MAX_PATH, L".\\setting.ini");
	wstring ws(value);
	changeObjName = string(ws.begin(), ws.end()); //wchar_t to std::string
	GetPrivateProfileString(L"Model", L"ShoeNameLeft", L"data/yaya_shoe2_left.obj", value, MAX_PATH, L".\\setting.ini");
	ws = wstring(value);
	changeObjNameL = string(ws.begin(), ws.end());
	GetPrivateProfileString(L"Model", L"FootName", L"data/new rice2.obj", value, MAX_PATH, L".\\setting.ini");
	ws = wstring(value);
	changeFootName = string(ws.begin(), ws.end());
	GetPrivateProfileString(L"Model", L"FootNameLeft", L"data/new rice2_left.obj", value, MAX_PATH, L".\\setting.ini");
	ws = wstring(value);
	changeFootNameL = string(ws.begin(), ws.end());

	// [Calibration]
	GetPrivateProfileString(L"Calibration", L"calibX", L"0", value, MAX_PATH, L".\\setting.ini");
	ws = wstring(value);
	calibX = atof(string(ws.begin(), ws.end()).c_str());
	GetPrivateProfileString(L"Calibration", L"calibY", L"0", value, MAX_PATH, L".\\setting.ini");
	ws = wstring(value);
	calibY = atof(string(ws.begin(), ws.end()).c_str());
	GetPrivateProfileString(L"Calibration", L"calibZ", L"0", value, MAX_PATH, L".\\setting.ini");
	ws = wstring(value);
	calibZ = atof(string(ws.begin(), ws.end()).c_str());
	GetPrivateProfileString(L"Calibration", L"calibF", L"0", value, MAX_PATH, L".\\setting.ini");
	ws = wstring(value);
	calibF = atof(string(ws.begin(), ws.end()).c_str());

	// [Perspective]
	GetPrivateProfileString(L"Perspective", L"fovy", L"0", value, MAX_PATH, L".\\setting.ini");
	ws = wstring(value);
	fovy = atof(string(ws.begin(), ws.end()).c_str());
	GetPrivateProfileString(L"Perspective", L"aspect", L"0", value, MAX_PATH, L".\\setting.ini");
	ws = wstring(value);
	aspect = atof(string(ws.begin(), ws.end()).c_str());

	// [Hint]
	GetPrivateProfileString(L"Hint", L"hintLX", L"0", value, MAX_PATH, L".\\setting.ini");
	ws = wstring(value);
	hintLX = atoi(string(ws.begin(), ws.end()).c_str());
	GetPrivateProfileString(L"Hint", L"hintLY", L"0", value, MAX_PATH, L".\\setting.ini");
	ws = wstring(value);
	hintLY = atoi(string(ws.begin(), ws.end()).c_str());
	GetPrivateProfileString(L"Hint", L"hintRX", L"0", value, MAX_PATH, L".\\setting.ini");
	ws = wstring(value);
	hintRX = atoi(string(ws.begin(), ws.end()).c_str());
	GetPrivateProfileString(L"Hint", L"hintRY", L"0", value, MAX_PATH, L".\\setting.ini");
	ws = wstring(value);
	hintRY = atoi(string(ws.begin(), ws.end()).c_str());

	// [Floor]
	GetPrivateProfileString(L"Floor", L"VectorX", L"0", value, MAX_PATH, L".\\setting.ini");
	ws = wstring(value);
	floor_Vector.x = atof(string(ws.begin(), ws.end()).c_str());
	GetPrivateProfileString(L"Floor", L"VectorY", L"0", value, MAX_PATH, L".\\setting.ini");
	ws = wstring(value);
	floor_Vector.y = atof(string(ws.begin(), ws.end()).c_str());
	GetPrivateProfileString(L"Floor", L"VectorZ", L"0", value, MAX_PATH, L".\\setting.ini");
	ws = wstring(value);
	floor_Vector.z = atof(string(ws.begin(), ws.end()).c_str());
	GetPrivateProfileString(L"Floor", L"VectorW", L"0", value, MAX_PATH, L".\\setting.ini");
	ws = wstring(value);
	floor_Vector.w = atof(string(ws.begin(), ws.end()).c_str());
	
}

void WriteINIFile()
{
	// [Calibration]
	wstringstream wss;
	wss << calibX;
	WritePrivateProfileString(L"Calibration", L"calibX ", wss.str().c_str(), L".\\setting.ini");

	wss.str(L"");
	wss.clear();
	wss << calibY;
	WritePrivateProfileString(L"Calibration", L"calibY", wss.str().c_str(), L".\\setting.ini");

	wss.str(L"");
	wss.clear();
	wss << calibZ;
	WritePrivateProfileString(L"Calibration", L"calibZ", wss.str().c_str(), L".\\setting.ini");

	wss.str(L"");
	wss.clear();
	wss << calibF;
	WritePrivateProfileString(L"Calibration", L"calibF", wss.str().c_str(), L".\\setting.ini");

	// [Perspective]
	wss.str(L"");
	wss.clear();
	wss << fovy;
	WritePrivateProfileString(L"Perspective", L"fovy", wss.str().c_str(), L".\\setting.ini");

	wss.str(L"");
	wss.clear();
	wss << aspect;
	WritePrivateProfileString(L"Perspective", L"aspect", wss.str().c_str(), L".\\setting.ini");

	// [Hint]
	wss.str(L"");
	wss.clear();
	wss << hintLX;
	WritePrivateProfileString(L"Hint", L"hintLX", wss.str().c_str(), L".\\setting.ini");

	wss.str(L"");
	wss.clear();
	wss << hintLY;
	WritePrivateProfileString(L"Hint", L"hintLY", wss.str().c_str(), L".\\setting.ini");

	wss.str(L"");
	wss.clear();
	wss << hintRX;
	WritePrivateProfileString(L"Hint", L"hintRX", wss.str().c_str(), L".\\setting.ini");

	wss.str(L"");
	wss.clear();
	wss << hintRY;
	WritePrivateProfileString(L"Hint", L"hintRY", wss.str().c_str(), L".\\setting.ini");

	// [Floor]
	wss.str(L"");
	wss.clear();
	wss << floor_Vector.x;
	WritePrivateProfileString(L"Floor", L"VectorX", wss.str().c_str(), L".\\setting.ini");

	wss.str(L"");
	wss.clear();
	wss << floor_Vector.y;
	WritePrivateProfileString(L"Floor", L"VectorY", wss.str().c_str(), L".\\setting.ini");

	wss.str(L"");
	wss.clear();
	wss << floor_Vector.z;
	WritePrivateProfileString(L"Floor", L"VectorZ", wss.str().c_str(), L".\\setting.ini");

	wss.str(L"");
	wss.clear();
	wss << floor_Vector.w;
	WritePrivateProfileString(L"Floor", L"VectorW", wss.str().c_str(), L".\\setting.ini");

	//WritePrivateProfileString(L"Calibration", L"calibF", wss.str().c_str(), L".\\setting.ini");
}

void adjustModelScale(int model, int size){
	if (model == 0){
		foot->Scale(1 / modelScale, 1 / modelScale, 1 / modelScale);
		footL->Scale(1 / modelScale, 1 / modelScale, 1 / modelScale);
	}
	else if (model == 1){
		//obj->Scale(1 / modelScale, 1 / modelScale, 1 / modelScale);
		hide->Scale(1 / modelScale, 1 / modelScale, 1 / modelScale);
		//objL->Scale(1 / modelScale, 1 / modelScale, 1 / modelScale);
		hideL->Scale(1 / modelScale, 1 / modelScale, 1 / modelScale);
	}

	if (size == 1)
		modelScale += .01;
	else
		modelScale -= .01;

	if (model == 0){
		foot->Scale(modelScale, modelScale, modelScale);
		footL->Scale(modelScale, modelScale, modelScale);
	}
	else if (model == 1){
		//obj->Scale(modelScale, modelScale, modelScale);
		hide->Scale(modelScale, modelScale, modelScale);
		//objL->Scale(modelScale, modelScale, modelScale);
		hideL->Scale(modelScale, modelScale, modelScale);
	}
}
void adjustModelTranslate(int dx, int dy, int dz){
	obj->Translate(dx, dy, dz);
	objL->Translate(dx, -dy, dz);
}

int main(int argc, char** argv){	

	thread socketThread(skt);
	init();
	ReadINIFile();
	wchar_t szDir[256];
	GetCurrentDirectory(256, szDir);
	wstring ws(szDir);
	// your new String
	string cDir(ws.begin(), ws.end());
	// Show String
	//cout << cDir << endl;	

	//---------------------------------------

	/*string sfoot = cDir + "\\data\\new rice2.obj";
	//string sfoot = cDir + "\\data/Foot/Foot-42-R-Low-Test.obj";
	string sobj = cDir + "\\data/Shoes/TF-8-R.obj";
	string shide = cDir + "\\data\\Remove\\Remove-R.obj";
	string sfootL = cDir + "\\data\\new rice2_left.obj";
	string sobjL = cDir + "\\data/Shoes/TF-8-L.obj";
	string shideL = cDir + "\\data\\Remove\\Remove-L.obj";*/

	/*
	//char footname[] = "data/Shoes/nike3.obj"; //純鞋子 3000面 砍後面
	//char footname[] = "data/Shoes/nikeR+foot2_low.obj"; // 鞋+足 3000面 沒砍後面 X
	//char footname[] = "data/Shoes/nikeR+foot2_low1.obj"; // 鞋+足 3000面 砍後面
	char footname[] = "data/Shoes/nikeR_low.obj"; // 鞋 3000面 砍後面(多)
	
	//char footname[] = "data/new rice2.obj"; // 原始足模
	//char objname[] = "data/Shoes/nike1.obj"; // 純鞋子
	//char objname[] = "data/Shoes/logo2_tx.obj"; // logo+texture
	char objname[] = "data/Shoes/Sole.obj"; // sole
	//char objname[] = "data/Shoes/nikeR_low_hide.obj";
	char hidename[] = "data/Shoes/nikeR_low_hide.obj";
	//char hidename[] = "data/Remove/Remove-R.obj";*/

	// foot model
	char footname[] = "data/new rice2.obj"; // 原始足模
	char footnameL[] = "data/new rice2_left.obj";

	// occlusion model
	//char hidename[] = "data/Remove/Remove-R-revise-specialtop.obj";
	char hidename[] = "data/Remove/Remove-R-sb.obj";
	//char hidenameL[] = "data/Remove/Remove-L-revise-specialtop.obj";
	char hidenameL[] = "data/Remove/Remove-L-sb.obj";

	// shoe model
	//char objname[] = "data/Shoes/TF-5-R.obj";
	//char objnameL[] = "data/Shoes/TF-5-L.obj";

	char objname[] = "data/Shoes/TF-4-R.obj";
	char objnameL[] = "data/Shoes/TF-4-L.obj";

	//char objnameL[] = "data/puma115105_left.obj";

	// memory locate
	img = (color3uc *)malloc(1920 * 1080 * sizeof(color3uc));
	d_img = (color3uc *)malloc(1920 * 1080 * sizeof(color3uc));
	depth = (point3f*)malloc(1920 * 1080 * sizeof(point3f));
	imgTest = (color3uc *)malloc(1388 * 1080 * sizeof(color3uc));

	// load obj. model
	obj = new Mesh();
	obj->ReadOBJ(objname);
	foot = new Mesh();
	foot->ReadOBJ(footname);
	hide = new Mesh();
	hide->ReadOBJ(hidename);
	tempObj = new Mesh();
	objL = new Mesh();
	objL->ReadOBJ(objnameL);
	footL = new Mesh();
	footL->ReadOBJ(footnameL);
	hideL = new Mesh();
	hideL->ReadOBJ(hidenameL);

	// scale
	float footScale = 0.95, shoeScale = 0.85;
/*	foot->Scale(footScale+0.05*1, footScale+0.05*1, footScale+0.05*1);
	obj->Scale(shoeScale+0.13*1, shoeScale+0.13*1, shoeScale+0.13*1);
	hide->Scale(footScale-0.01*1, footScale-0.01*1, footScale-0.01*1);*/
	
	obj->Scale(shoeScale, shoeScale, shoeScale);
	obj->Translate(-15.0, -0.0, -5.0);

	objL->Scale(shoeScale, shoeScale, shoeScale);
	objL->Translate(-15.0, -0.0, -5.0);
	
	foot->Scale(footScale, footScale, footScale);
	footL->Scale(footScale, footScale, footScale);
	
	hide->Scale(footScale, footScale*1.1, footScale);
	hide->Translate(-40.0, -28.0, -20.0);
	
	hideL->Scale(footScale, footScale*1.1, footScale);
	hideL->Translate(-35.0, 0.0, -10.0);

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
	GL_Setup();
	//glutInitWindowSize(1920, 1080);
	
	int startingX, startingY, windowHeight, windowWidth;
	if (argc == 2){
		startingY = atoi(argv[1]);
		glutInitWindowPosition(0, startingY);
		glutInitWindowSize(2160, 1215);
	}else if(argc == 3){
		startingY = atoi(argv[1]);
		startingX = atoi(argv[2]);
		glutInitWindowPosition(startingX, startingY);
		glutInitWindowSize(2160, 1215);
	}
	else if (argc == 5){
		startingY = atoi(argv[1]);
		startingX = atoi(argv[2]);
		glutInitWindowPosition(startingX, startingY);
		windowHeight = atoi(argv[3]);
		windowWidth = atoi(argv[4]);
		glutInitWindowSize(windowWidth, windowHeight);
	}
	else{
		//glutInitWindowPosition(0, 2160);
		//glutInitWindowSize(2160, 1920);
		//glutInitWindowPosition(0, 0);
		glutInitWindowSize(1920, 1080);
	}
	
	WinNumber = glutCreateWindow("Virtual Try-on");
	HWND win_handle = FindWindow(0, L"Virtual Try-on");
	SetWindowLong(win_handle, GWL_STYLE, (GetWindowLong(win_handle, GWL_STYLE) | WS_MAXIMIZE));
	ShowWindowAsync(win_handle, SW_SHOWMAXIMIZED);
	//SetWindowStyle();
	//glutFullScreen(); //解析度才會是1920*1080

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_LIGHTING);

	// default camera
	
	glutDisplayFunc(display);
	glutIdleFunc(idle);
	glutKeyboardFunc(keyboard);
	glutMouseFunc(Mouse);
	glutTimerFunc(1000, &timer, 0);
	//glutReshapeFunc(reshape);

	glutMainLoop();
	return 0;

}
void init(){
	// 1. Sensor related code
	cout << "Try to get default sensor" << endl;
	{
		if (GetDefaultKinectSensor(&pSensor) != S_OK)
		{
			cerr << "Get Sensor failed" << endl;
			return ;
		}

		cout << "Try to open sensor" << endl;
		if (pSensor->Open() != S_OK)
		{
			cerr << "Can't open sensor" << endl;
			return ;
		}
	}

	// 2. Color related code
	cout << "Try to get color source" << endl;
	{
		// Get frame source
		IColorFrameSource* pFrameSource = nullptr;
		if (pSensor->get_ColorFrameSource(&pFrameSource) != S_OK)
		{
			cerr << "Can't get color frame source" << endl;
			return ;
		}

		// Get frame description
		cout << "get color frame description\n" << endl;
		IFrameDescription* pFrameDescription = nullptr;
		if (pFrameSource->get_FrameDescription(&pFrameDescription) == S_OK)
		{
			pFrameDescription->get_Width(&iColorWidth);
			pFrameDescription->get_Height(&iColorHeight);

			uColorPointNum = iColorWidth * iColorHeight;
			uColorBufferSize = uColorPointNum * 4 * sizeof(BYTE);

			pCSPoints = new CameraSpacePoint[uColorPointNum];
			pColorBuffer = new BYTE[4 * uColorPointNum];
			printf("uColorPointNum = %d\n", uColorPointNum);
		}
		pFrameDescription->Release();
		pFrameDescription = nullptr;

		// get frame reader
		cout << "Try to get color frame reader" << endl;
		if (pFrameSource->OpenReader(&pColorFrameReader) != S_OK)
		{
			cerr << "Can't get color frame reader" << endl;
			return ;
		}

		// release Frame source
		cout << "Release frame source" << endl;
		pFrameSource->Release();
		pFrameSource = nullptr;
	}

	// 3. Depth related code
	cout << "Try to get depth source" << endl;
	{
		// Get frame source
		IDepthFrameSource* pFrameSource = nullptr;
		if (pSensor->get_DepthFrameSource(&pFrameSource) != S_OK)
		{
			cerr << "Can't get depth frame source" << endl;
			return ;
		}

		// Get frame description
		cout << "get depth frame description\n" << endl;
		IFrameDescription* pFrameDescription = nullptr;
		if (pFrameSource->get_FrameDescription(&pFrameDescription) == S_OK)
		{
			int	iDepthWidth = 0,
				iDepthHeight = 0;
			pFrameDescription->get_Width(&iDepthWidth);
			pFrameDescription->get_Height(&iDepthHeight);
			uDepthPointNum = iDepthWidth * iDepthHeight;
			pDepthBuffer = new UINT16[uDepthPointNum];
			printf("uDepthPointNum = %d\n", uDepthPointNum);
		}
		pFrameDescription->Release();
		pFrameDescription = nullptr;

		// get frame reader
		cout << "Try to get depth frame reader" << endl;
		if (pFrameSource->OpenReader(&pDepthFrameReader) != S_OK)
		{
			cerr << "Can't get depth frame reader" << endl;
			return ;
		}
		// release Frame source
		cout << "Release frame source" << endl;
		pFrameSource->Release();
		pFrameSource = nullptr;
	}

	// 4. Coordinate Mapper
	if (pSensor->get_CoordinateMapper(&pCoordinateMapper) != S_OK)
	{
		cerr << "get_CoordinateMapper failed" << endl;
		return ;
	}	

	// 5. Body related code
	{
		//IBodyFrameReader*	pBodyFrameReader = nullptr;
		// Get frame source
		IBodyFrameSource* pFrameSource = nullptr;
		if (pSensor->get_BodyFrameSource(&pFrameSource) != S_OK)
		{
			cerr << "Can't get depth frame source" << endl;
			return;
		}

		// get frame reader
		//cout << "Try to get depth frame reader" << endl;
		if (pFrameSource->OpenReader(&pBodyFrameReader) != S_OK)
		{
			cerr << "Can't get depth frame reader" << endl;
			return;
		}
		// release Frame source
		//cout << "Release frame source" << endl;
		pFrameSource->Release();
		pFrameSource = nullptr;
	}
}

// glut idle function
time_t totalTime = 0;
void idle()
{
	//time_t tick,tick2,captureTime;
	// Get floor data
	//captureTime = clock();
	IBodyFrame* pBodyFrame = NULL;
	if (pBodyFrameReader->AcquireLatestFrame(&pBodyFrame) == S_OK)
	{
		//INT64 nTime = 0;
		//pBodyFrame->get_RelativeTime(&nTime);
		//pBodyFrame->get_FloorClipPlane(&floor_vector);
		//printf("f_x = %f, f_y = %f, f_z = %f, f_w = %f\n", floor_vector.x, floor_vector.y, floor_vector.z, floor_vector.w);
		pBodyFrame->Release();
		pBodyFrame = nullptr;
	}	
	//cout << "Get floor data: " << clock() - tick << endl;
	// Read color data
	//tick = clock();
	IColorFrame* pCFrame = nullptr;
	if (pColorFrameReader->AcquireLatestFrame(&pCFrame) == S_OK)
	{
		pCFrame->CopyConvertedFrameDataToArray(uColorBufferSize, pColorBuffer, ColorImageFormat_Rgba);
		
		pCFrame->Release();
		pCFrame = nullptr;
	}
	//std::cout << "Read color data: " << (float)(clock() - tick) / CLOCKS_PER_SEC << endl;
	// Read depth data
	//tick = clock();
	IDepthFrame* pDFrame = nullptr;
	if (pDepthFrameReader->AcquireLatestFrame(&pDFrame) == S_OK)
	{
		pDFrame->CopyFrameDataToArray(uDepthPointNum, pDepthBuffer);

		pDFrame->Release();
		pDFrame = nullptr;

		// map to camera space
		pCoordinateMapper->MapColorFrameToCameraSpace(uDepthPointNum, pDepthBuffer, uColorPointNum, pCSPoints);
		
		
		//tick2 = clock();
		//ByteArraryToColorImg(); //0.01 sec
		//AlignDepthMaps(); //0.04 sec =>0.03(thread)
		thread mThread(ByteArraryToColorImg);
		thread mThread2(AlignDepthMaps);
		mThread.join();
		mThread2.join();
		//std::cout << (float)(clock() - tick2) / CLOCKS_PER_SEC << endl;
			
			
	}
	//initail or draw
	if (isinit == 0){
		isinit = 1;
		InitTracker(depth, 1920, 1080);
		for (int i = 0; i < window_height; i++){ // 1080 down
			for (int j = 0; j < window_width; j++){ // 1920 right
				//index1 = (1079 - i) * window_width + j;
				depth[(1079 - i) * window_width + j].x = 0;
				depth[(1079 - i) * window_width + j].y = 0;
				depth[(1079 - i) * window_width + j].z = 0;
				zValue[(1079 - i) * window_width + j] = 0;
			}
		}
	}
	else{		
		// redraw
		FPS++;
		glutPostRedisplay();
		//std::cout << "Total time: " << (float)(clock() - totalTime) / CLOCKS_PER_SEC << " sec" << endl;
		//totalTime = clock();
	}
	//std::cout << "Read depth data: " << clock() - tick << endl;
	//std::cout << "Kinect capture time: " << (float)(clock() - captureTime) / CLOCKS_PER_SEC << " sec" << endl;
}

void AlignDepth(int threadIdx)
{

	floor_vector.x = floor_Vector.x;
	floor_vector.y = floor_Vector.y;
	floor_vector.z = floor_Vector.z;
	floor_vector.w = floor_Vector.w;

	//int ltY = 400, ltX = 660, ltH = 500, ltW = 600;
	float x, y, z;
	int index, index2;
	//校正參數
	calibX += distancex;
	calibY += distancey;
	calibZ += distancez;
	calibF += distancef;
	distancex = 0;
	distancey = 0;
	distancez = 0;
	distancef = 0;
	/*Vector4 floor_Vector;
	floor_Vector.x = 1.79118942e-005;
	floor_Vector.y = 0.895879447;
	floor_Vector.z = -0.444297343;
	floor_Vector.w = 0.666652203;*/

	for (int i = 1080 / alignThreadNum*threadIdx; i < 1080 / alignThreadNum*(threadIdx + 1); i++){ // 1080 down
		for (int j = 0; j < 1920; j++){ // 1920 right
			//if (isinit == 0 || (isinit != 0 && (_i >= ltY && _i < ltY + ltH && _j >= ltX && _j < ltX + ltW)))
			{
				index = i * 1920 + j;
				index2 = (1079 - i) * 1920 + j;
				x = pCSPoints[index].X + calibX;
				y = pCSPoints[index].Y + calibY;
				z = pCSPoints[index].Z + calibZ;

				//if (z <= 0 || (x * floor_Vector.x + y * floor_Vector.y + z * floor_Vector.z + floor_Vector.w < calibF))//if (z <= 0 ) // z <= 0 是不正常的深度點  後面判斷式是將地板深度點拿掉
				if (z <= 0 || (x * floor_vector.x + y * floor_vector.y + z * floor_vector.z + floor_vector.w < calibF))//if (z <= 0 ) // z <= 0 是不正常的深度點  後面判斷式是將地板深度點拿掉
				{
					depth[index2].x = 0;
					depth[index2].y = 0;
					depth[index2].z = 0;
					zValue[index2] = 0;
				}
				else
				{
					depth[index2].x = x * 1000;
					depth[index2].y = y * 1000;
					depth[index2].z = -z * 1000;
					zValue[index2] = z * 1000;
				}
			}
		}
	}
}
void AlignDepthMaps()
{
	for (int i = 0; i < alignThreadNum; i++)
	{
		alignThread[i] = thread(AlignDepth, i);
	}
	for (int i = 0; i < alignThreadNum; i++)
	{
		alignThread[i].join();
	}
}
void ByteArraryToColorImg()
{
	for (int i = 0; i < window_height; i++){ // 1080 down
		for (int j = 0; j < window_width; j++){ // 1920 right
			//memcpy(img, pColorBuffer, sizeof(BYTE) * 3 * window_height*window_width);
			img[(1079 - i) * window_width + j].r = pColorBuffer[4 * (i * window_width + j)];
			img[(1079 - i) * window_width + j].g = pColorBuffer[4 * (i * window_width + j) + 1];
			img[(1079 - i) * window_width + j].b = pColorBuffer[4 * (i * window_width + j) + 2];
		}
	}
}
void ChangeShoes()
{
	obj->ReadOBJ(changeObjName);
}
void ChangeShoesL()
{
	objL->ReadOBJ(changeObjNameL);
}
void display()
{
	
	//time_t displayTime;
	//displayTime = clock();
	glViewport(0, 0, 1920, 1080);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//gluPerspective(52.0f + viewval, 1.285 + ratioval, 0.01f, 20000.0f);
	gluPerspective(fovy, aspect, 0.05f, 20000.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
	if (testmode == 1){
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		DrawBackgroundUC(img);	

		//-----------Tracking----------------
		thread threadR, threadL;
		switch (mode){
		case 1:
			SetTrackerTexture(img);
			SetTrackerDepth(zValue);
			right_feet.TrackColor(FT_GREEN, 20, m);
			break;
		case 2:			
			threadR = thread(RightFeetTracking);
			threadL = thread(LeftFeetTracking);
			threadR.join();
			threadL.join();
			draw = 2;
			break;
		case 3:
			SetTrackerDepth(zValue);
			right_feet.test(m);
			break;
		case 4:			
			selectpoints.n = 0;
			selectpoints.pt = (Point3f *)malloc(6 * sizeof(Point3f));
			selectpointsL.n = 0;
			selectpointsL.pt = (Point3f *)malloc(6 * sizeof(Point3f));
			mode = 5;
			break;
		case 5:
			if (selectpoints.n >= 6)
			{		
				right_feet.TrackSelect(selectpoints, m);
				mode = 6; 
			}else{	
				

				glPointSize(10.0);
				glBegin(GL_POINTS);

				glColor3f(1.0, 0.0, 0.0);				
				for (int i = 0; i < selectpoints.n; i++)
					glVertex3f(selectpoints.pt[i].x, selectpoints.pt[i].y, selectpoints.pt[i].z);

				glEnd();
				glPointSize(1.0);
			}
			break;
		case 6:
			RightFeetTracking();
			draw = 2;
			if (selectpointsL.n >= 6)
			{
				
			//	right_feet.TrackSelect(selectpoints, m);
		/*		m[0] = 0.0; m[1] = 1.0; m[2] = 0.0;
				m[4] = 0.0; m[5] = 0.0; m[6] = 1.0;
				m[8] = 1.0; m[9] = 0.0; m[10] = 0.0;
				m[3] = selectpoints.pt[0].x;
				m[7] = selectpoints.pt[0].y;
				m[11] = selectpoints.pt[0].z;*/
				left_feet.TrackSelect(selectpointsL, mL);
			/*	for (int i = 0; i < 16; i++){
					if (isnan(m[i]) || isinf(m[i]) || isnan(mL[i]) || isinf(mL[i]))
					{
						cout << "error: nan" << endl;
						selectpoints.n = 0;
						selectpointsL.n = 0;
						mode = 4;
						break;
					}
				}*/
				/*cout << m[0] << "," << m[1] << "," << m[2] << "," << m[3] << endl;
				cout << m[4] << "," << m[5] << "," << m[6] << "," << m[7] << endl;
				cout << m[8] << "," << m[9] << "," << m[10] << "," << m[11] << endl;
				cout << m[12] << "," << m[13] << "," << m[14] << "," << m[15] << endl;*/
				//cout << m[3] << "," << m[7] << "," << m[11] << "," << m[15] << endl;
				
				mode = 2; //mode = 2;
			}
			else{
		/*		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
				glBegin(GL_POLYGON);
				glColor3f(1.0, 0.0, 0.0);
				glVertex3f(-50, -150, -350);
				glVertex3f(-50, -50, -350);
				glVertex3f(-150, -50, -350);
				glVertex3f(-150, -150, -350);
				glEnd();
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);*/
				glPointSize(10.0);
				glBegin(GL_POINTS);

				glColor3f(1.0, 0.0, 0.0);
				for (int i = 0; i < selectpointsL.n; i++)
					glVertex3f(selectpointsL.pt[i].x, selectpointsL.pt[i].y, selectpointsL.pt[i].z);

				glEnd();
				glPointSize(1.0);
			}
			break;
		case 7:
			
		/*	glPointSize(30.0);
			glBegin(GL_POINTS);

			glColor3f(1.0, 1.0, 0.0);
			glVertex3f(depth[(1079 - 680) * 1920 + 1200].x, depth[(1079 - 680) * 1920 + 1200].y, depth[(1079 - 680) * 1920 + 1200].z);
			glEnd();*/

			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glLineWidth(5.0);
			//glBegin(GL_POLYGON);
			if (depth[(1079 - hintLY) * 1920 + hintLX].z != 0)
				glColor3f(0.05*(initTimer), 1, 0.05*(initTimer));				
			else
				glColor3f(1.0, 0.0, 0.0);
			{
				glPointSize(15.0);
				glBegin(GL_POINTS);
				glVertex3f(depth[(1079 - hintLY) * 1920 + hintLX].x, depth[(1079 - hintLY) * 1920 + hintLX].y, depth[(1079 - hintLY) * 1920 + hintLX].z);
				glEnd();
			}
			/*glVertex3f(40 + distancex, -80 + distancey, -350 + distancez);
			glVertex3f(40 + distancex, -10 + distancey, -350 + distancez);
			glVertex3f(110 + distancex, -10 + distancey, -350 + distancez);
			glVertex3f(110 + distancex, -80 + distancey, -350 + distancez);*/
			//glEnd();
			//glBegin(GL_POLYGON);
			if (depth[(1079 - hintRY) * 1920 + hintRX].z != 0)
				glColor3f(0.05*(initTimer), 1, 0.05*(initTimer));
			else
				glColor3f(1.0, 0.0, 0.0);
			{
				glPointSize(15.0);
				glBegin(GL_POINTS);
				glVertex3f(depth[(1079 - hintRY) * 1920 + hintRX].x, depth[(1079 - hintRY) * 1920 + hintRX].y, depth[(1079 - hintRY) * 1920 + hintRX].z);
				glEnd();
			}
			/*glVertex3f(30, -80, -350);
			glVertex3f(30, -10, -350);
			glVertex3f(-40, -10, -350);
			glVertex3f(-40, -80, -350);
			glEnd();*/
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			if (depth[(1079 - hintLY) * 1920 + hintLX].z != 0 && depth[(1079 - hintRY) * 1920 + hintRX].z != 0)
				initTimer++;
			else
				initTimer = 0;
			if (initTimer > 20)
			{
				right_feet.TrackInitSolution(m, depth[(1079 - hintRY) * 1920 + hintRX]);
				left_feet.TrackInitSolution(mL, depth[(1079 - hintLY) * 1920 + hintLX]);				
				mode = 2;
			}
			glPointSize(1.0);
			break;

		case 8:
			// position
			/*RT_Paramf theta[1];
			float m1[3][4];

			theta[0].x = -12;
			theta[0].y = -137;
			theta[0].z = -1100;
			theta[0].a = 0.568555;
			theta[0].b = -0.743519;
			theta[0].c = -2.108131;

			CompRTMatf(m1, theta[0]);
			memcpy(m, m1[0], 12 * sizeof(float));

			draw = 2;*/
			break;
		default:			
			break;
		}


		//--------------------- set up occlusion improvement -------------------------
		// define foot center point
		float footmodelpoint[3] = { /*-69.9890*/ -140, -13.2558, /*35.4902*/ 35 };
		float footmodelpointL[3] = { /*-69.9890*/ -140, 13.2558, /*35.4902*/ 35 };

		// left foot
		float footmodelpointTopL2[3] = { -175, 35, 120 };
		float footmodelpointTopL[3] = { -155, 13.2558, 70 };
		float footanglepointL[3] = { -127.5, 13.2558, 10 };

		// right foot
		float FootModel_Top_R[3] = { -155, -13.2558, 70 };
		float FootAnglePoint_R[3] = { -127.5, -13.2558, 10 };

		// trans
		float Newfootmodelpoint[3] = { m[0] * footmodelpoint[0] + m[1] * footmodelpoint[1] + m[2] * footmodelpoint[2] + m[3],
			m[4] * footmodelpoint[0] + m[5] * footmodelpoint[1] + m[6] * footmodelpoint[2] + m[7],
			m[8] * footmodelpoint[0] + m[9] * footmodelpoint[1] + m[10] * footmodelpoint[2] + m[11] };
		float NewfootmodelpointL[3] = { mL[0] * footmodelpointL[0] + mL[1] * footmodelpointL[1] + mL[2] * footmodelpointL[2] + mL[3],
			mL[4] * footmodelpointL[0] + mL[5] * footmodelpointL[1] + mL[6] * footmodelpointL[2] + mL[7],
			mL[8] * footmodelpointL[0] + mL[9] * footmodelpointL[1] + mL[10] * footmodelpointL[2] + mL[11] };
		float NewfootmodelpointTopL2[3] = { mL[0] * footmodelpointTopL2[0] + mL[1] * footmodelpointTopL2[1] + mL[2] * footmodelpointTopL2[2] + mL[3],
			mL[4] * footmodelpointTopL2[0] + mL[5] * footmodelpointTopL2[1] + mL[6] * footmodelpointTopL2[2] + mL[7],
			mL[8] * footmodelpointTopL2[0] + mL[9] * footmodelpointTopL2[1] + mL[10] * footmodelpointTopL2[2] + mL[11] };


		float NewfootmodelpointTopL[3] = { mL[0] * footmodelpointTopL[0] + mL[1] * footmodelpointTopL[1] + mL[2] * footmodelpointTopL[2] + mL[3],
			mL[4] * footmodelpointTopL[0] + mL[5] * footmodelpointTopL[1] + mL[6] * footmodelpointTopL[2] + mL[7],
			mL[8] * footmodelpointTopL[0] + mL[9] * footmodelpointTopL[1] + mL[10] * footmodelpointTopL[2] + mL[11] };
		float NewfootanglepointL[3] = { mL[0] * footanglepointL[0] + mL[1] * footanglepointL[1] + mL[2] * footanglepointL[2] + mL[3],
			mL[4] * footanglepointL[0] + mL[5] * footanglepointL[1] + mL[6] * footanglepointL[2] + mL[7],
			mL[8] * footanglepointL[0] + mL[9] * footanglepointL[1] + mL[10] * footanglepointL[2] + mL[11] };

		float Trans_FootModel_Top_R[3] = { m[0] * FootModel_Top_R[0] + m[1] * FootModel_Top_R[1] + m[2] * FootModel_Top_R[2] + m[3],
			m[4] * FootModel_Top_R[0] + m[5] * FootModel_Top_R[1] + m[6] * FootModel_Top_R[2] + m[7],
			m[8] * FootModel_Top_R[0] + m[9] * FootModel_Top_R[1] + m[10] * FootModel_Top_R[2] + m[11] };
		float Trans_FootAnglePoint_R[3] = { m[0] * FootAnglePoint_R[0] + m[1] * FootAnglePoint_R[1] + m[2] * FootAnglePoint_R[2] + m[3],
			m[4] * FootAnglePoint_R[0] + m[5] * FootAnglePoint_R[1] + m[6] * FootAnglePoint_R[2] + m[7],
			m[8] * FootAnglePoint_R[0] + m[9] * FootAnglePoint_R[1] + m[10] * FootAnglePoint_R[2] + m[11] };

		// ---------- Draw foot occlusion point ----------

		footanglepoint = DepthOfAnkle(NewfootmodelpointTopL[0], NewfootmodelpointTopL[1], NewfootmodelpointTopL[2], 120, NewfootanglepointL[0], NewfootanglepointL[1], NewfootanglepointL[2], 
			0.0, 0.057f, -0.009f, 0.06f);
		footanglepoint2 = DepthOfAnkle(Trans_FootModel_Top_R[0], Trans_FootModel_Top_R[1], Trans_FootModel_Top_R[2], 120, Trans_FootAnglePoint_R[0], Trans_FootAnglePoint_R[1], Trans_FootAnglePoint_R[2], 
			1.0, 0.057f, -0.009f, 0.06f);

		
		//-----------Drawing model----------------
		Light();
		glEnable(GL_LIGHTING);
		CompGLMatf(glm, m);
		glPushMatrix();
		glMultMatrixf(glm);
		
		if (draw != 5)
		{
			glEnable(GL_BLEND);
			glBlendFunc(GL_ZERO, GL_ONE);
			//hide->Draw(MESH_FLAT);
			glBlendFunc(GL_ONE, GL_ZERO);
		}
		switch (draw){
		case 1:
			obj->Draw(MESH_VERTEX | MESH_MATERIAL | MESH_NORMAL | MESH_TEXTURE);
			break;
		case 2:			
			if (shoeNum == 0){
				obj->Draw(MESH_SMOOTH | MESH_MATERIAL | MESH_NORMAL | MESH_TEXTURE);
			}
			else{
				foot->Draw(MESH_SMOOTH | MESH_MATERIAL | MESH_NORMAL | MESH_TEXTURE);
			}
			break;
		case 3:
			obj->Draw(MESH_WIRE | MESH_SMOOTH | MESH_MATERIAL | MESH_NORMAL | MESH_TEXTURE);
			break;
		case 5:
			foot->Draw(MESH_WIRE | MESH_SMOOTH | MESH_NORMAL);
			break;
		case 6:			
			break;
		case 4:
			break;
		}	
		glPopMatrix();
		CompGLMatf(glmL, mL);
		glPushMatrix();
		glMultMatrixf(glmL);

		if (draw != 5)
		{
			glEnable(GL_BLEND);
			glBlendFunc(GL_ZERO, GL_ONE);
			//hideL->Draw(MESH_FLAT);
			glBlendFunc(GL_ONE, GL_ZERO);
		}
		switch (draw){
		case 1:
			objL->Draw(MESH_VERTEX | MESH_MATERIAL | MESH_NORMAL | MESH_TEXTURE);
			break;
		case 2:
			if (shoeNum == 0){
				objL->Draw(MESH_SMOOTH | MESH_MATERIAL | MESH_NORMAL | MESH_TEXTURE);
			}
			else{
				footL->Draw(MESH_SMOOTH | MESH_MATERIAL | MESH_NORMAL | MESH_TEXTURE);
			}
			break;
		case 3:
			objL->Draw(MESH_WIRE | MESH_SMOOTH | MESH_MATERIAL | MESH_NORMAL | MESH_TEXTURE);
			break;
		case 5:
			footL->Draw(MESH_WIRE | MESH_SMOOTH | MESH_NORMAL);
			break;
		case 6:
			break;
		case 4:
			break;
		}
		glPopMatrix();
		glDisable(GL_LIGHTING);

		glutSwapBuffers();
	}
	else if (testmode == 2){
		// clear previous screen
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		//-----繪圖(深度圖)-----
		DrawDepthPointcloud(img, depth);


		// swap buffer
		glutSwapBuffers();
	}

	if (needChangeShoe)
	{
		obj->Clear();
		objL->Clear();
		ExitTracker(); //free some memory to load jpeg texture
		/*time_t tc;
		tc = clock();*/
		thread changeThread1(ChangeShoes);
		thread changeThread2(ChangeShoesL);
		changeThread1.join();
		changeThread2.join();
		//std::cout << "ChangeShoes: " << (float)(clock() - tc) / CLOCKS_PER_SEC << std::endl;
		isinit = 0; //re-initial tracker
		needChangeShoe = false;
		
	}

	//std::cout << "Display time: " << (float)(clock() - displayTime) / CLOCKS_PER_SEC << " sec" << endl;
}
void RightFeetTracking()
{
	//SetTrackerTexture(img);
	//time_t t;
	//t = clock();
	SetTrackerDepth(zValue);
	//std::cout << "SetTrackerDepth: " << clock() - t << std::endl;
	mbuf.v = foot->GetVertex(0);
	mbuf.vn = foot->GetNormal(0);
	mbuf.f = (struct _face_buffer *)foot->GetTriangle(0);
	mbuf.n = foot->GetVertexSize();
	mbuf.nf = foot->GetTriangleSize();
	//t = clock();
	right_feet.TrackICP(mbuf, m, icpe, FT_RIGHT);
	//std::cout << "TrackICP: " << clock() - t << std::endl;
}
void LeftFeetTracking()
{
	//SetTrackerTexture(img);
	//time_t t;
	//t = clock();
	SetTrackerDepth(zValue);
	//std::cout << "SetTrackerDepth: " << clock() - t << std::endl;
	mbufL.v = footL->GetVertex(0);
	mbufL.vn = footL->GetNormal(0);
	mbufL.f = (struct _face_buffer *)footL->GetTriangle(0);
	mbufL.n = footL->GetVertexSize();
	mbufL.nf = footL->GetTriangleSize();
	//t = clock();
	left_feet.TrackICP(mbufL, mL, icpeL, FT_LEFT);
	//std::cout << "TrackICP: " << clock() - t << std::endl;
}
void timer(int value)
{
	glutTimerFunc(1000, &timer, 0);
	printf("FPS = %d\n", FPS);
	FPS = 0;	
}
void Mouse(int button, int state, int x, int y)
{
	if (state)
	{
		printf("%d %d		%f %f %f \n", x, y, \
			depth[(1079 - y) * 1920 + x].x, \
			depth[(1079 - y) * 1920 + x].y, \
			depth[(1079 - y) * 1920 + x].z);

		if (depth[(1079 - y) * 1920 + x].z != 0)
		{
			if (mode == 5)
			{
				selectpoints.pt[selectpoints.n].x = depth[(1079 - y) * 1920 + x].x;
				selectpoints.pt[selectpoints.n].y = depth[(1079 - y) * 1920 + x].y;
				selectpoints.pt[selectpoints.n].z = depth[(1079 - y) * 1920 + x].z;
				selectpoints.n++;
			}
			else if (mode == 6)
			{
				selectpointsL.pt[selectpointsL.n].x = depth[(1079 - y) * 1920 + x].x;
				selectpointsL.pt[selectpointsL.n].y = depth[(1079 - y) * 1920 + x].y;
				selectpointsL.pt[selectpointsL.n].z = depth[(1079 - y) * 1920 + x].z;
				selectpointsL.n++;
			}
		}
	}

}
void reshape(int w, int h)
{
	//glViewport(0, 0, 640, 480);
	glViewport(0, 0, w, h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45.0, (float)w / h, 0.05, 50.0);
	//gluPerspective(45.0f, 1.333f, 0.1f, 20000.0f);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	
}
void keyboard(unsigned char key, int x, int y)
{
	switch (key)
	{	
	case '1':
		testmode = 1;
		break;
	case '2':
		testmode = 2;		
		break;
	case '3':
		mode = 4;
		draw = 6;
		break;
	case '4':
		shoeNum = 1 - shoeNum;
		break;
	case '5':
		draw = 5;
		break;
	case '6':
		draw = 2;
		break;
	case '7':
		// save calibration parameters
		WriteINIFile();
		break;
	case '0':
		mode = 7;
		initTimer = 0;
		m[11] = -100.0;
		mL[11] = -100.0;
		break;
	case 'w':	
		switch (adjustmode){
		case 1:
			hintLY += 10;
			break;
		case 2:
			hintRY += 10;
			break;
		case 3:
			adjustModelTranslate(0, 0, 5);
			break;
		default:
			distancey += .001;
			break;
		}
		/*if (adjustmode == 1){
			hintLY += 10;
		}
		else if (adjustmode == 2){
			hintRY += 10;
		}
		else{
			distancey += .001;
		}*/
		break;
	case 's':	
		switch (adjustmode){
		case 1:
			hintLY -= 10;
			break;
		case 2:
			hintRY -= 10;
			break;
		case 3:
			adjustModelTranslate(0, 0, -5);
			break;
		default:
			distancey -= .001;
			break;
		}
		/*if (adjustmode == 1){
			hintLY -= 10;
		}
		else if (adjustmode == 2){
			hintRY -= 10;
		}
		else{
			distancey -= .001;
		}*/
		break;
	case 'a':	
		switch (adjustmode){
		case 1:
			hintLX -= 10;
			break;
		case 2:
			hintRX -= 10;
			break;
		case 3:
			adjustModelTranslate(0, -5, 0);
			break;
		case 4:
			adjustModelScale(0, -1);
			break;
		case 5:
			adjustModelScale(1, -1);
			break;
		default:
			distancex -= .001;
			break;
		}
		/*if (adjustmode == 1){
			hintLX -= 10;
		}
		else if (adjustmode == 2){
			hintRX -= 10;
		}
		else{
			distancex -= .001;
		}*/
		break;
	case 'd':		
		switch (adjustmode){
		case 1:
			hintLX += 10;
			break;
		case 2:
			hintRX += 10;
			break;
		case 3:
			adjustModelTranslate(0, 5, 0);
			break;
		case 4:
			adjustModelScale(0, 1);
			break;
		case 5:
			adjustModelScale(1, 1);
			break;
		default:
			distancex += .001;
			break;
		}
		/*if (adjustmode == 1){
			hintLX += 10;
		}
		else if (adjustmode == 2){
			hintRX += 10;
		}
		else{
			distancex += .001;
		}*/
		break;
	case 'z':
		if (adjustmode == 3)
			adjustModelTranslate(-5, 0, 0);
		else
			distancez -= .001;
		break;
	case 'x':
		if (adjustmode == 3)
			adjustModelTranslate(5, 0, 0);
		else
			distancez += .001;
		break;
	case 'r':
		distancef += .001;
		break;
	case 'f':
		distancef -= .001;
		break;
	case 'p':
		adjustmode = (++adjustmode) % 6;
		if (adjustmode == 1){
			printf("=====Adjust Left Foot Hint Position=====\n");
		}
		else if (adjustmode == 2){
			printf("=====Adjust Right Foot Hint Position=====\n");
		}
		else if (adjustmode == 3){
			printf("=====Adjust Shoe Position=====\n");
		}
		else if (adjustmode == 4){
			printf("=====Adjust Foot Size=====\n");
		}
		else if (adjustmode == 5){
			printf("=====Adjust Shoe Size=====\n");
		}
		break;
	case '[':
		adjustModelScale(0, 1);
		break;
	case ']':
		adjustModelScale(0, -1);
		break;
	case ';':
		adjustModelScale(1, 1);
		break;
	case '\'':
		adjustModelScale(1, -1);
		break;
	case 't':
		fovy += 0.1;
		//viewval += 1;
		break;
	case 'g':
		fovy -= 0.1;
		//viewval -= 1;
		break;
	case 'y':
		aspect += .005;
		//ratioval += .005;
		break;
	case 'h':
		aspect -= .005;
		//ratioval -= .005;
		break;
	case 27:
		glDisable(GL_LIGHT0);
		glDisable(GL_LIGHTING);
		glDisable(GL_DEPTH_TEST);
		glutDestroyWindow(WinNumber);
		exit(0);
		break;	
	default:
		printf("%d", key);
		break;
	}
	//printf("distancex = %f, distancey = %f, distancez = %f, viewval = %f, ratioval = %f\n", distancex, distancey, distancez, viewval, ratioval);
	//printf("Floor_vector : x = %f, y = %f, z = %f, w = %f\n", floor_vector.x, floor_vector.y, floor_vector.z, floor_vector.w);
	printf("rate = %f, modelScale = %f\n", rate, modelScale);
}
