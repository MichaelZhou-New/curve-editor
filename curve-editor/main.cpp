#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <limits>
#include <functional>

#include "Geometry.h"
#include "GLDebug.h"
#include "Log.h"
#include "ShaderProgram.h"
#include "Shader.h"
#include "Texture.h"
#include "Window.h"

#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"



// By Dylan Leclair (UCID: 30068119) with some external resources listed as needed.

// We gave this code in one of the tutorials, so leaving it here too
void updateGPUGeometry(GPU_Geometry &gpuGeom, CPU_Geometry &cpuGeom, int selindex) {
	if (selindex >= 0) {
		cpuGeom.cols[selindex] = glm::vec3(0.5f, 0.0f, 0.5f);
	}

	gpuGeom.bind();
	gpuGeom.setVerts(cpuGeom.verts);
	gpuGeom.setCols(cpuGeom.cols);



}


void updateGPUGeometry(GPU_Geometry& gpuGeom, CPU_Geometry const& cpuGeom) {
	gpuGeom.bind();
	gpuGeom.setVerts(cpuGeom.verts);
	gpuGeom.setCols(cpuGeom.cols);


}

// EXAMPLE CALLBACKS
class Assignment3 : public CallbackInterface {

public:
	Assignment3(int screenWidth, int screenHeight) : screenDimensions(screenWidth, screenHeight)
	{
	}


	void doTheViewPipeline(ShaderProgram& shader, float deltaTime) {

		// getting the location of V and P from shader 
		// (i did not go thru the pain of using a model)
		GLint VLoc = glGetUniformLocation(shader, "V");
		GLint PLoc = glGetUniformLocation(shader, "P");

		// initialize P and V so they can be set accordingly 
		glm::mat4 P, V;


		if (is3d) {

			// See cursor position callback for more on this stuff

			const float cameraSpeed = 1.0f * deltaTime;

			if (isWDown)
				cameraPos += cameraSpeed * cameraFront;
			if (isSDown)
				cameraPos -= cameraSpeed * cameraFront;
			if (isADown)
				cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
			if (isDDown)
				cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;


			// calculation of view and perspective matrices for 3D viewer

			V = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
			P = glm::perspective(glm::radians(65.0f), 1.0f, 0.1f, 100.0f);

		}
		else {

			// update 2d params
			if (isWDown) { // zoom 
				
				if (zoom - 0.05 > 0.0f) { // make sure they dont go past the curve...
					zoom -= 0.05;
				}
				
			}
			if (isSDown) { // zoom out
				zoom += 0.05;
			}

			// calculation of view and perspective matrices for curve editor (orthogonal seemed to fit 2D better)

			V = glm::lookAt(glm::vec3(0.0f, 0.0f, 1.0f), glm::vec3(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
			P = glm::ortho(-1 * zoom, zoom, -1 * zoom, zoom, 0.1f, 100.0f);



		}

		// Upload V and P to the shader to avoid catastophe

		glUniformMatrix4fv(VLoc, 1, false, glm::value_ptr(V));
		glUniformMatrix4fv(PLoc, 1, false, glm::value_ptr(P));




	}

	virtual void keyCallback(int key, int scancode, int action, int mods) {


		// TOGGLES

		// Switches between curve editor and 3d viewer when K is pressed
		if (key == GLFW_KEY_K && action == GLFW_PRESS) {
			if (is3d) {
				is3d = false;
				// resets camera to original position to make things nice and simple
				cameraPos = glm::vec3(0.0f, 0.0f, 3.0f);
				cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
				cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
			}
			else {
				is3d = true;
			}
		}


		// Toggles between bezier and b-spline when B is pressed
		if (key == GLFW_KEY_B && action == GLFW_PRESS) {
			if (bspline) {
				bspline = false;
			}
			else {
				bspline = true;
			}
		}


		// Toggles wireframe when J is pressed
		if (key == GLFW_KEY_J && action == GLFW_PRESS) {
			
			if (is3d) {
				if (wireframe) {
					wireframe = false;
				}
				else {
					wireframe = true;
				}
			}
			

		}

		// Toggles surface control points
		if (key == GLFW_KEY_P && action == GLFW_PRESS) {
			if (showSurfacePoints) {
				showSurfacePoints = false;
			}
			else {
				showSurfacePoints = true;
			}
		}

		// Curve point selection

		if (key == GLFW_KEY_TAB && action == GLFW_PRESS) {
			tab = true;
		}
		if (key == GLFW_KEY_TAB && action == GLFW_RELEASE) {
			tab = false;
		}

		// Curve clear key

		if (key == GLFW_KEY_C && action == GLFW_PRESS) {
			clear = true;
		}



		// CAMERA / ZOOM CONTROL HANDLING
		// If a button is pressed, set the appropriate bool to true until it is released.

		if (key == GLFW_KEY_W && action == GLFW_PRESS) {
			isWDown = true;
		}
		if (key == GLFW_KEY_S && action == GLFW_PRESS) {
			isSDown = true;
		}

		if (key == GLFW_KEY_A && action == GLFW_PRESS) {
			isADown = true;
		}
		if (key == GLFW_KEY_D && action == GLFW_PRESS) {
			isDDown = true;
		}

		if (key == GLFW_KEY_W && action == GLFW_RELEASE) {
			isWDown = false;
		}
		if (key == GLFW_KEY_S && action == GLFW_RELEASE) {
			isSDown = false;
		}

		if (key == GLFW_KEY_A && action == GLFW_RELEASE) {
			isADown = false;
		}
		if (key == GLFW_KEY_D && action == GLFW_RELEASE) {
			isDDown = false;
		}


		// 3D scene selection
		// I realize now I have done the if statement to check for this once
		// on the outside instead but I've spent too much time on this
		// to refactor it too.

		if (key == GLFW_KEY_1 && action == GLFW_PRESS) {

			if (is3d) { 
				scene3d = 1;
			}
			
		}
		if (key == GLFW_KEY_2 && action == GLFW_PRESS) {
			if (is3d) {
				scene3d = 2;
			}
		}
		if (key == GLFW_KEY_3 && action == GLFW_PRESS) {
			if (is3d) {
				scene3d = 3;
			}
		}


	}


	virtual void mouseButtonCallback(int button, int action, int mods) {


		if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
			// when the left mouse is clicked, alter it's bool to true
			leftMousePressed = true;
		}


		if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
			// on click of right mouse button
			rightMousePressed = true;
		}




	}
	virtual void cursorPosCallback(double xpos, double ypos) {

		// I got a lot of this (everything but formatting / organization) from
		// https://learnopengl.com/Getting-started/Camera
		// I figured since we never really talked about cameras and such in class
		// that it wasn't a super important thing for us to slave over.

		// This was partly motivated by my whopping C+ in math 211
		// (I don't know how I've made it this far either)

		// GIVES SCREEN COORDINATES. MUST BE CONVERTED TO WORLD COORDINATES
		// X AND Y START IN TOP LEFT (see mouseGL)
		xScreenPos = xpos;
		yScreenPos = ypos;


		if (is3d) {

			if (first)
			{
				lastX = xpos;
				lastY = ypos;

				//std::cout << "can't catch me" << std::endl;
				first = false;
				//cameraFront = glm::vec3(0.0f, 0.0f, -1.0f)
			}
			else {
				float xoffset = xpos - lastX;
				float yoffset = lastY - ypos;
				lastX = xpos;
				lastY = ypos;

				float sensitivity = 0.1f;

				xoffset *= sensitivity;
				yoffset *= sensitivity;

				yaw += xoffset;
				pitch += yoffset;

				if (pitch > 89.0f)
					pitch = 89.0f;
				if (pitch < -89.0f)
					pitch = -89.0f;

				glm::vec3 direction;
				direction.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
				direction.y = sin(glm::radians(pitch));
				direction.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
				cameraFront = glm::normalize(direction);
			}

		}

	}

	// I considered using this but decided it was too cursed
	virtual void scrollCallback(double xoffset, double yoffset) {
	}

	// please don't change the window size I don't have the heart to see how broken everything is if it's changed.
	virtual void windowSizeCallback(int width, int height) {
		// The CallbackInterface::windowSizeCallback will call glViewport for us
		CallbackInterface::windowSizeCallback(width,  height);
	}

	// performs necessary transformations to get coordinates from cursor position.
	// this is from tutorial!
	glm::vec2 mouseGL() {
		glm::vec2 startingVec(xScreenPos, yScreenPos);
		glm::vec2 shiftedVec = startingVec + glm::vec2(0.5f);

		glm::vec2 scaledToZeroOne = shiftedVec / glm::vec2(screenDimensions); // component wise division (x divided by x, y divided by y, ...)

		glm::vec2 flippedY = glm::vec2(scaledToZeroOne.x, 1.0f - scaledToZeroOne.y);

		glm::vec2 final = flippedY * 2.0f - glm::vec2(1.0f);

		final.x = final.x * zoom;
		final.y = final.y * zoom;

		return final;
	}

	// A nice long list of getters / setters for all my parameters.

	bool isLeftMouseDown() {
		return leftMousePressed;
	}
	bool isRightMouseDown() {
		return rightMousePressed;
	}

	bool isTabDown() {
		return tab;
	}

	void refreshMouse() {
		leftMousePressed = false;
		rightMousePressed = false;
	}


	bool isBSpline() {
		return bspline;
	}

	bool isView3D() {
		return is3d;
	}
	bool isWireframe() {
		return wireframe;
	}

	bool isClear() {
		return clear;
	}

	void refreshClear() {
		clear = false;
	}

	int get3dScene() {
		return scene3d;
	}

	bool isShowSurfacePoints() {
		return showSurfacePoints;
	}

private:
	bool leftMousePressed = false; // was the left mouse button pressed
	bool rightMousePressed = false; // was the right mouse button pressed
	glm::vec2 screenDimensions; // self explanatory
	double xScreenPos = 0;
	double yScreenPos = 0;
	bool tab = false; // is tab down?

	float zoom = 1.2; // default zoom (for curve editor)
	bool is3d = false; // is the user viewing 3d object or editing curve?
	bool bspline = false; // use bspline / chaikin subdivision if true

	// bool params for camera
	bool isWDown = false;
	bool isSDown = false;
	bool isADown = false;
	bool isDDown = false;
	bool first = true; 

	// if true, render the control points for the surface being viewed.
	bool showSurfacePoints = false;

	// encodes the scene (1 = shape created from curve, 2 = first surface, 3 = second surface)
	int scene3d = 1;

	bool wireframe = false; // whether or not to use wireframe

	// more camera parameters...
	glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 3.0f);
	glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
	glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
	float yaw = -90.0f;
	float pitch = 0.0f;
	float lastX = 400.0f;
	float lastY = 400.0f;

	// whether or not the curve should be cleared 
	bool clear = false;
};





// returns the deCasteljou Q(u) (some point during the bezier curve), u is between 0 and 1
glm::vec3 deCastel(const std::vector<glm::vec3>& control, int d, float u) {

	std::vector<glm::vec3> nums = control;

	for (int i = 1; i <= d; ++i) {
		for (int j = 0; j <= (d - i); ++j) {
			nums[j] = (1 - u) * nums[j] + u * nums[j + 1];
			
		}
	}

	return nums[0];
}

// Generates vertices and colours for a deCasteljou curve into CPU geometry
CPU_Geometry generateCasteljou(const CPU_Geometry& cpoints) {

	CPU_Geometry curve;

	for (float i = 0; i <= 1; i += 0.01f) {
		// generate the point using deCastel
		glm::vec3 point = deCastel(cpoints.verts, cpoints.verts.size() - 1, i);
		//std::cout << point << std::endl;
		curve.verts.push_back(point);
		//curve.cols.push_back(glm::vec3(0.0f, 1.0f, 1.0f)); // black
	}

	// get the last point in there?
	curve.verts.push_back(cpoints.verts[cpoints.verts.size() - 1]);
	

	return curve;
}

// Generates vertices and colours for a deCasteljou curve into CPU geometry
std::vector<glm::vec3> generateCasteljouVec(const CPU_Geometry& cpoints) {

	std::vector<glm::vec3> curve;

	for (float i = 0; i <= 1; i += 0.05f) {
		// generate the point using deCastel
		glm::vec3 point = deCastel(cpoints.verts, cpoints.verts.size() - 1, i);
		//std::cout << point << std::endl;
		curve.push_back(point); // add it to our collection
		//curve.cols.push_back(glm::vec3(0.0f, 1.0f, 1.0f)); // black
	}

	// get the last point in there?
	curve.push_back(cpoints.verts[cpoints.verts.size() - 1]);


	return curve;
}


// Generates vertices  for a deCasteljou curve into a vector
// Returns vector of points on curve (from which geometry must be generated)
std::vector<glm::vec3> generateCasteljouVec(const std::vector<glm::vec3>& cpoints) {

	std::vector<glm::vec3> curve;

	for (float i = 0; i <= 1; i += 0.05f) {
		// generate the point using deCastel
		glm::vec3 point = deCastel(cpoints, cpoints.size() - 1, i);
		//std::cout << point << std::endl;
		curve.push_back(point);
		//curve.cols.push_back(glm::vec3(0.0f, 1.0f, 1.0f)); // black
	}

	// get the last point in there?
	curve.push_back(cpoints[cpoints.size() - 1]);


	return curve;
}

// Uses chaikin subdivision to generate B-Spline curve
std::vector<glm::vec3> chaikin(std::vector<glm::vec3> control, int n) {
	std::vector<glm::vec3> result;


	result.push_back((control[0]));
	result.push_back(((0.5f * control[0]) + (0.5f * control[1])));

	// periodic mask in middle

	for (int i = 1; i < n - 2; i++) {
		result.push_back(0.75f * control[i] + 0.25f * control[i + 1]);
		result.push_back(0.25f * control[i] + 0.75f * control[i + 1]);
	}

	// mask at end

	result.push_back(0.5f * control[n - 2] + 0.5f *control[n-1]);
	result.push_back(control[n-1]);

	return result;
}

// Uses chaikin subdivision to generate B-Spline curve
// Driver for chaikin
// Returns vector of points on curve (from which geometry must be generated)
std::vector<glm::vec3> generateBSpline(std::vector<glm::vec3> control, int iterations) {


	for (int i = 0; i < iterations; ++i) {
		control = chaikin(control, control.size());
	}

	return control;
	
}

// Used as a helper for revolution - finds the radius from the Y axis of a point
float radiusY(glm::vec3 point) {

	return sqrt(point.x * point.x + point.z * point.z);

}

// generates the revolution surface of from the control points of a given curve
// loop increments can be modified to change quality of object
// Returns vector of points on curve (from which geometry must be generated)
std::vector < std::vector< glm::vec3 >> revolveY(const std::vector<glm::vec3>& control) {

	std::vector < std::vector< glm::vec3 >> object;

	for (float u = 0; u < 360; u += 5) {
		float rads = glm::radians(u);

		std::vector<glm::vec3> rev;
		for (int v = 0; v < control.size(); v+=2) {

			float r = radiusY(control[v]);

			rev.push_back(glm::vec3(r * cos(rads), control[v].y, r * sin(rads)));
		}
		object.push_back(rev);

	}

	return object;

}

// Creates the triangles needed to properly render the object, from the points representing it.
// Returns the appropriate geometry which should be rendered with GL_TRIANGLES
CPU_Geometry generateRevolutionGeom(std::vector < std::vector< glm::vec3 >> revArray) {
	CPU_Geometry doughnut;

	for (int i = 0; i < revArray.size() - 1; ++i) {
		std::vector<glm::vec3> current = revArray[i];
		std::vector<glm::vec3> next = revArray[i + 1];

		int elemsize = current.size() - 1;

		for (int j = 0; j < elemsize; j++) {
			// draw a face
			doughnut.verts.push_back(current[j]);
			doughnut.verts.push_back(next[j]);
			doughnut.verts.push_back(current[j + 1]);

			doughnut.verts.push_back(current[j + 1]);
			doughnut.verts.push_back(next[j + 1]);
			doughnut.verts.push_back(next[j]);

		}

	}
	// WRAP BEGGINING AND END

	int daddy = revArray.size() - 1;

	auto current = revArray[daddy];
	auto next = revArray[0];

	int elemsize = current.size() - 1;

	for (int j = 0; j < elemsize; j++) {
		// draw a face
		doughnut.verts.push_back(current[j]);
		doughnut.verts.push_back(next[j]);
		doughnut.verts.push_back(current[j + 1]);

		doughnut.verts.push_back(current[j + 1]);
		doughnut.verts.push_back(next[j + 1]);
		doughnut.verts.push_back(next[j]);

	}

	// Reset the colors to green
	doughnut.cols.clear();
	doughnut.cols.resize(doughnut.verts.size(), glm::vec3{ 1.0, 0.0, 1.0 });

	return doughnut;
}

// Finds the euclidean distance between two points.
float euclidDistance(glm::vec3 a, glm::vec3 b) {
	float dist;

	float dx = a.x - b.x;
	float dy = a.y - b.y;

	dist = pow(dx, 2) + pow(dy, 2);
	dist = sqrt(dist);

	return dist;
}

// Given a 2d array of control points, generates the 2d array of finer points from which to generate the surface.
// uses the deCasteljou algorithm to do subdivision // smoothing
std::vector<std::vector<glm::vec3>> calculateSurfaceCasteljou(const std::vector<std::vector<glm::vec3>>& control) {

	// based off of the psuedocode in lecture

	std::vector<std::vector<glm::vec3>> R;

	// calculate each R row and add this to result vector, return it
	for (int rowIndex = 0; rowIndex < control.size(); rowIndex++) {

		std::vector<glm::vec3> r = generateCasteljouVec(control[rowIndex]);
		R.push_back(r);

	}

	std::vector<std::vector<glm::vec3>> Q;


	for (int columnIndex = 0; columnIndex < R[0].size(); columnIndex++) {

		// accumulate the column values

		std::vector<glm::vec3> column;

		for (int row = 0; row < R.size(); row++) {

			// go down each COLUMN in R

			column.push_back(R[row][columnIndex]);


		}

		// do casteljau on the column

		std::vector<glm::vec3> q = generateCasteljouVec(column);

		// add this to Q (field of final points representing surface)
		Q.push_back(q);
	}

	return Q; // must be stitched together to create a quad


}


// Given a 2d array of control points, generates the 2d array of finer points from which to generate the surface.
// uses the chaikin algorithm to do subdivision // smoothing (b-spline curves)
std::vector<std::vector<glm::vec3>> calculateSurfaceChaikin(const std::vector<std::vector<glm::vec3>>& control) {


	std::vector<std::vector<glm::vec3>> R;


	// calculate each R row and add this to result vector, return it
	for (int rowIndex = 0; rowIndex < control.size(); rowIndex++) {

		std::vector<glm::vec3> r = generateBSpline(control[rowIndex], 4);
		R.push_back(r);

	}

	std::vector<std::vector<glm::vec3>> Q;


	for (int columnIndex = 0; columnIndex < R[0].size(); columnIndex++) {

		// accumulate the column values

		std::vector<glm::vec3> column;

		for (int row = 0; row < R.size(); row++) {


			// Add each column value to column...

			column.push_back(R[row][columnIndex]);


		}

		// do casteljau on the column

		std::vector<glm::vec3> q = generateBSpline(column, 4);

		// add this to Q (field of final points representing surface)
		Q.push_back(q);
	}

	return Q; // must be stitched together to create a quad


}

// I used this to generate random y-values for a flat surface to generate 
// the tensor flow surface control points.
// from https://stackoverflow.com/questions/13408990/how-to-generate-random-float-number-in-c
float float_rand(float min, float max)
{
	float scale = rand() / (float)RAND_MAX; /* [0, 1.0] */
	return min + scale * (max - min);      /* [min, max] */
}

// Generates control points for a surface
// The idea behind this was to create a flat surface on X-Z plane
// then change the elevation of points at random to create 
// random surfaces every iteration
// i wasn't about to hard code that ;-; 
// and i get way cooler surfaces...
std::vector<std::vector<glm::vec3>> generateSurfaceControlPoints(int numRows, int numCols) {

	std::vector<std::vector<glm::vec3>> control;

	glm::vec3 start = glm::vec3(1.0f, 0.0f, 1.0f);

	for (int i = 0; i < numRows; i++) {
		float xIncrement = -1 * i / static_cast<float>(numRows);

		// for each point
		glm::vec3 start = glm::vec3(1.0f + xIncrement, 0.0f, 1.0f);

		std::vector<glm::vec3> row;
		for (float j = 0; j < numCols; j++) {

			glm::vec3 point = start;

			float zIncrement = -1 * j / static_cast<float>(numCols);

			point.z += zIncrement;

			// if an edge, use y = 0
			if (j == 0 || j == numCols - 1) {
				point.y = 0.0f;
			}
			else { // otherwise, use a random y!
				point.y = float_rand(-0.6f, 0.6f);
			}


			row.push_back(point);

		}

		// collect the rows
		control.push_back(row);
		

	}

	return control;

}


// Stiches a 2d array of points together to create the geometry needed
// to draw the shape with GL_TRIANGLES.
CPU_Geometry generateStitchedGeom(std::vector < std::vector< glm::vec3 >> revArray) {
	CPU_Geometry doughnut; // i tried to do this with a circle at first to draw a doughnut
	// i was also hungry... 
	// the name stuck.


	// basically, for each row....

	for (int i = 0; i < revArray.size() - 1; ++i) {
		std::vector<glm::vec3> current = revArray[i];
		std::vector<glm::vec3> next = revArray[i + 1];

		int elemsize = current.size() - 1;

		// we take the next row and draw a bunch of squares connecting the two rows.
		// but computers dont like squares so we actually drew a triangle instead
		// each vertex in the loop below is some corner of the triangle.
		// I had a diagram I drew but I can't really comment that in

		for (int j = 0; j < elemsize; j++) {
			// draw a face
			doughnut.verts.push_back(current[j]);
			doughnut.verts.push_back(next[j]);
			doughnut.verts.push_back(current[j + 1]);

			doughnut.verts.push_back(current[j + 1]);
			doughnut.verts.push_back(next[j + 1]);
			doughnut.verts.push_back(next[j]);

		}

	}

	// Reset the colors to green
	doughnut.cols.clear();
	doughnut.cols.resize(doughnut.verts.size(), glm::vec3{ 1.0, 0.0, 1.0 });

	return doughnut;
}


// Collapses a 2D array of control points into a 1D array so that it may be rendered
// idk if this was necessary but like 
// it works
CPU_Geometry collectSurfaceControlPoints(const std::vector < std::vector< glm::vec3 >>& control) {
	CPU_Geometry cpts;

	for (auto row = control.begin(); row < control.end(); row++) {
		std::vector< glm::vec3 > r = *row;
		for (auto col = r.begin(); col < r.end(); col++) {
			glm::vec3 point = *col;
			cpts.verts.push_back(point);
		}

	}

	cpts.cols.clear();
	cpts.cols.resize(cpts.verts.size(), glm::vec3{ 1.0, 0.0, 0.0 });
	return cpts;
}

int main() {
	Log::debug("Starting main");

	// WINDOW
	glfwInit();
	Window window(800, 800, "Curve Editor"); // can set callbacks at construction if desired


	GLDebug::enable();

	// CALLBACKS
	auto a3 = std::make_shared<Assignment3>(800.0f,800.0f);
	window.setCallbacks(a3);


	ShaderProgram shader("shaders/test.vert", "shaders/test.frag");

	// The current CPU_Geometry and GPU_Geometry classes are defined in
	// Geometry.h/Geometry.cpp They will work for this assignment, but for some of
	// the bonuses you may have to modify them.

	// initial points for curve cuz why not
	CPU_Geometry square; // in hindsight this was a bad name
	square.verts.push_back(glm::vec3{-0.5, 0.5, 0});
	square.verts.push_back(glm::vec3{-0.5, -0.5, 0});
	square.verts.push_back(glm::vec3{0.5, -0.5, 0});
	square.verts.push_back(glm::vec3{0.5, 0.5, 0});

	// SQUARE actually represents the vector of control points from which to create a curve.

	// generating all the GPU geometry... (from boilerplate)

	square.cols.clear();
	square.cols.resize(square.verts.size(), glm::vec3{1.0, 0.0, 0.0});
	GPU_Geometry pointsGPUGeom;
	updateGPUGeometry(pointsGPUGeom, square);

	// Reset the colors to green
	square.cols.clear();
	square.cols.resize(square.verts.size(), glm::vec3{0.0, 1.0, 0.0});

	GPU_Geometry linesGPUGeom;
	updateGPUGeometry(linesGPUGeom, square);


	// the actual curve created from the control points of square.
	CPU_Geometry curve;
	curve.verts = generateBSpline(square.verts, 12);
	curve.cols.resize(curve.verts.size(), glm::vec3(0.0f, 0.0f, 0.0f));
	GPU_Geometry bcurve;
	updateGPUGeometry(bcurve, curve);



	// SURFACES & STUFF
	// you can kind of tell I knew what I was doing for this part because it's not just spaghetti
	// I did this all at once because we didn't need to alter the surfaces (thank god)

	srand(time(nullptr));

	// Generate control points for two somewhat-random surfaces
	std::vector < std::vector< glm::vec3 >> surface1control = generateSurfaceControlPoints(4, 4);
	std::vector < std::vector< glm::vec3 >> surface2control = generateSurfaceControlPoints(4, 4);

	// the actual finer surface points for them go into these bad boys
	std::vector < std::vector< glm::vec3 >> surface1, surface1spline, surface2, surface2spline;

	// create the geometry (triangles) for both, using both bezier and bspline 
	CPU_Geometry surf1geom, surf1bgeom, surf2geom, surf2bgeom;

	// create the bezier surface for each
	surface1 = calculateSurfaceCasteljou(surface1control);
	surface2 = calculateSurfaceCasteljou(surface2control);

	// Gather the control points for each surface into a 1D vector
	CPU_Geometry surface1Control = collectSurfaceControlPoints(surface1control);
	CPU_Geometry surface2Control = collectSurfaceControlPoints(surface2control);

	// Create the GPU_Geometry to draw the control points
	GPU_Geometry surfPoints1, surfPoints2;
	updateGPUGeometry(surfPoints1, surface1Control);
	updateGPUGeometry(surfPoints2, surface2Control);

	// create the b-spline surface for each
	surface1spline = calculateSurfaceChaikin(surface1control);
	surface2spline = calculateSurfaceChaikin(surface2control);


	// Now that we have the surfaces, we need to create ALL the triangles
	surf1geom = generateStitchedGeom(surface1);
	surf1bgeom = generateStitchedGeom(surface1spline);

	surf2geom = generateStitchedGeom(surface2);
	surf2bgeom = generateStitchedGeom(surface2spline);


	// Upload these triangles to the GPU geometry
	GPU_Geometry surf1GPU, surf1bGPU, surf2GPU, surf2bGPU;
	updateGPUGeometry(surf1GPU, surf1geom);
	updateGPUGeometry(surf1bGPU, surf1bgeom);
	updateGPUGeometry(surf2GPU, surf2geom);
	updateGPUGeometry(surf2bGPU, surf2bgeom);

	// now all we have to do is choose which one of these to draw

	glPointSize(10.0f);

	// Note this call only work on some systems, unfortunately.
	// In order for them to work, you have to comment out line 60
	// If you're on a mac, you can't comment out line 60, so you
	// these will have no effect. :(
	// glLineWidth(5.0f);

	// RENDER LOOP

	// a bunch of parameters for input - directly mirror a parameter/function from callbacks.
	bool leftClicked = false;
	bool rightClicked = false;
	bool tabbed = false;
	bool lastTabFrame = false;
	bool isView3d = false;
	int selectionindex = -1;
	bool bspline;
	bool isWireframe = false;
	glm::vec2 mousePos = a3->mouseGL();
	glm::vec3 closestPoint = square.verts.at(0);
	float dist = euclidDistance(glm::vec3(mousePos, 0.0f), closestPoint);

	bool new3dGeom = true;
	bool isClear = false;
	bool showSurfacePoints;

	int scene = 1;

	std::vector < std::vector< glm::vec3 >> object; //2d array of 3d points (stores the array of points for curve-generated object)
	CPU_Geometry objectTriangles; // the cpu geometry for the curve-generated object (from stitching the points together)

	auto timeElapsed = glfwGetTime();
	while (!window.shouldClose()) {

		
		// update time stuff
		auto newTimeElapsed = glfwGetTime();
		auto dt = newTimeElapsed - timeElapsed;
		timeElapsed = newTimeElapsed;

		// clear the bools that must be cleared
		a3->refreshMouse();
		a3->refreshClear();
		
		// do i even need to...
		glfwPollEvents();

		// Carries over a bunch of the information from callbacks
		// I realize this could be designed better but :///
		leftClicked = a3->isLeftMouseDown();
		rightClicked = a3->isRightMouseDown();
		mousePos = a3->mouseGL();
		tabbed = a3->isTabDown();
		bspline = a3->isBSpline();
		isView3d = a3->isView3D();
		isWireframe = a3->isWireframe();
		isClear = a3->isClear();
		scene = a3->get3dScene();
		showSurfacePoints = a3->isShowSurfacePoints();


		// if we live and in 3d
		if (isView3d) {

			// create the 3d geometry from the curve, only if the curve has changed at all (or is being run for first time)
			if (new3dGeom) {
				object = revolveY(curve.verts);
				objectTriangles = generateRevolutionGeom(object);

				objectTriangles.cols.clear();
				objectTriangles.cols.resize(objectTriangles.verts.size(), glm::vec3{ 1.0, 0.0, 1.0 });
				new3dGeom = false;
			}

		}
		else {
			


			if (isClear) {
				square.verts.clear();
				square.verts.push_back(glm::vec3(0.0f, 0.0f, 0.0f));
			}
			else {
				// MOVE THE CLOSEST POINT
				if (tabbed) {

					new3dGeom = true; // if a vertex was selected, redraw 3d geometry when switching to it

					if (selectionindex == -1) { // no selection in storage, select a point

						dist = euclidDistance(glm::vec3(mousePos, 0.0f), square.verts[0]);

						// find the closest point
						for (int i = 0; i < square.verts.size(); i++) {
							float localdist = euclidDistance(glm::vec3(mousePos, 0.0f), square.verts[i]);
							if (localdist <= dist) {
								closestPoint = square.verts[i];
								dist = localdist;

								if (dist <= 0.1) // ensure that there is some sense of locality ;-; 
									selectionindex = i;
							}
						}
					}
					else { // otherwise, a point has already been selected
						// translate the point to follow the mouse

						square.verts[selectionindex] = glm::vec3(mousePos.x, mousePos.y, 0.0f);

						// draw the new curve

						if (bspline) {
							curve.verts = generateBSpline(square.verts, 4);
						}
						else {
							curve.verts = generateCasteljouVec(square);
						}

						curve.cols.clear();
						curve.cols.resize(curve.verts.size(), glm::vec3(0.0f, 0.0f, 0.0f));
						updateGPUGeometry(bcurve, curve);

						// Start with the lines
						square.cols.clear();
						square.cols.resize(square.verts.size(), glm::vec3{ 1.0, 0.0, 0.0 });
						updateGPUGeometry(pointsGPUGeom, square, selectionindex);

						// Reset the colors to green
						square.cols.clear();
						square.cols.resize(square.verts.size(), glm::vec3{ 0.0, 1.0, 0.0 });
						updateGPUGeometry(linesGPUGeom, square);


					}



				}
				else {
					selectionindex = -1; // nothing is being selected

					// adding points
					if (leftClicked) {
						glm::vec2 mousePos = a3->mouseGL();
						//std::cout << "x: " << mousePos.x << " y:" << mousePos.y << std::endl;
						glm::vec3 newPoint = glm::vec3(mousePos.x, mousePos.y, 0.0f);
						square.verts.push_back(newPoint);



					}

					if (rightClicked) {
						// want to look thru verts and remove the icky point

						// some initialization
						glm::vec2 mousePos = a3->mouseGL();
						glm::vec3 closestPoint = square.verts.at(0);
						float dist = euclidDistance(glm::vec3(mousePos, 0.0f), closestPoint);

						std::vector<glm::vec3>::iterator pos;

						// find the closest point
						for (auto i = square.verts.begin(); i < square.verts.end(); i++) {
							float localdist = euclidDistance(glm::vec3(mousePos, 0.0f), *i);
							if (localdist <= dist) {
								closestPoint = *i;
								dist = localdist;
								pos = i;
							}
						}

						// delete the closest point
						// perhaps add a condition here so that this only happens within a given range to avoid accidental / unintentional deletions?
						if (square.verts.size() > 1) { // lets not crash anything shall we
							square.verts.erase(pos);
						}




					}

					if (leftClicked || rightClicked) {
						new3dGeom = true; // makes it so that we draw the 3d geometry again if points added / deleted
					}


					// check which type of curve to draw
					if (bspline) {
						curve.verts = generateBSpline(square.verts, 4);
					}
					else {
						curve = generateCasteljou(square);
					}

					// update colours
					curve.cols.clear();
					curve.cols.resize(curve.verts.size(), glm::vec3(0.0f, 0.0f, 0.0f));
					updateGPUGeometry(bcurve, curve);

					// Start with the points
					square.cols.clear();
					square.cols.resize(square.verts.size(), glm::vec3{ 1.0, 0.0, 0.0 });
					updateGPUGeometry(pointsGPUGeom, square);


					// Then do the lines
					square.cols.clear();
					square.cols.resize(square.verts.size(), glm::vec3{ 0.0, 1.0, 0.0 });
					updateGPUGeometry(linesGPUGeom, square);

				}
			}

		}


		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_FRAMEBUFFER_SRGB);
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


		// wireframe op
		if (isWireframe) {
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		}
		else {
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		}

		shader.use();

		// performs the necessary view / projection transformations
		a3->doTheViewPipeline(shader, dt);

		// time to draw!

		if (isView3d) { // if we are in 3d mode

			int size = 0; // initializing size of shape to render

			GPU_Geometry triangles; // initializing geometry for shape from curve editor

			if (scene == 1) { // scene 1 = curve editor object
				
				updateGPUGeometry(triangles, objectTriangles);

				size = objectTriangles.verts.size();

				triangles.bind();
			}
			else if (scene == 2) { // scene 2 = first surface
				if (bspline) {
					surf1bGPU.bind();
					size = surf1bgeom.verts.size();
				}
				else {
					surf1GPU.bind();
					size = surf1geom.verts.size();
				}

			}
			else if (scene == 3) { // scene 3 = second surface
				if (bspline) {
					surf2bGPU.bind();
					size = surf1bgeom.verts.size();
				}
				else {
					surf2GPU.bind();
					size = surf1geom.verts.size();
				}
			}

			GLFWwindow* winRef = window.getWindow();
			glfwSetInputMode(winRef, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

			// draw the surface/object
			glDrawArrays(GL_TRIANGLES, 0, GLsizei(size));

			// if control points have been enabled for surface, draw them
			if (showSurfacePoints) {
				if (scene == 2) {

					//draw control points for first surface
					surfPoints1.bind();
					glDrawArrays(GL_POINTS, 0, GLsizei(surface1Control.verts.size()));

				}
				else if (scene == 3) {
					// draw control points for second surface
					surfPoints2.bind();
					glDrawArrays(GL_POINTS, 0, GLsizei(surface2Control.verts.size()));
				}
			}
			

		}
		else { // if in curve editor, draw curve, control points, connection between control points.
			GLFWwindow* winRef = window.getWindow(); 
			glfwSetInputMode(winRef, GLFW_CURSOR, GLFW_CURSOR_NORMAL);

			linesGPUGeom.bind();
			glDrawArrays(GL_LINE_STRIP, 0, GLsizei(square.verts.size()));

			pointsGPUGeom.bind();
			glDrawArrays(GL_POINTS, 0, GLsizei(square.verts.size()));

			bcurve.bind();
			glDrawArrays(GL_LINE_STRIP, 0, GLsizei(curve.verts.size()));

		}

		
		glDisable(GL_FRAMEBUFFER_SRGB); // disable sRGB for things like imgui

		window.swapBuffers();
	}

	glfwTerminate();
	return 0;
}
