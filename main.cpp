#include <glut.h>

void DisplayScene() {
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0, 0.0, 0.0);

	glBegin(GL_QUADS);
	glVertex2f(450.0, 400.0);
	glVertex2f(800.0, 400.0);
	glVertex2f(800.0, 600.0);
	glVertex2f(450.0, 600.0);
	glEnd();

	glBegin(GL_QUADS);
	glVertex2f(0.0, 400.0);
	glVertex2f(350.0, 400.0);
	glVertex2f(350.0, 600.0);
	glVertex2f(0.0, 600.0);
	glEnd();

	glBegin(GL_TRIANGLES);
	glVertex2f(0.0, 200.0);
	glVertex2f(200.0, 400.0);
	glVertex2f(0.0, 400.0);
	glEnd();

	glBegin(GL_TRIANGLES);
	glVertex2f(550.0, 400.0);
	glVertex2f(600.0, 250.0);
	glVertex2f(650.0, 400.0);
	glEnd();

	glFlush();
}

void myinit() {
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glColor3f(1.0, 0.0, 0.0);
	glPointSize(1.0);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, 800.0, 600.0, 0.0);
}

void main(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowSize(800, 600);
	glutInitWindowPosition(0, 0);
	glutCreateWindow("Points");
	glutDisplayFunc(DisplayScene);

	myinit();
	glutMainLoop();
}