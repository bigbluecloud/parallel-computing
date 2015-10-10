#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include <X11/Xlib.h> //X11 library headers
#include <X11/Xutil.h>
#include <X11/Xos.h>

#define X_RESN 1000 //X resolution
#define Y_RESN 1000 //Y resolution
#define REAL_MAX 2 //Real and imaginary ranges for the complex numbers
#define REAL_MIN -2
#define IMAG_MAX 2
#define IMAG_MIN -2

typedef struct { //Complex number struct
	float real, imaginary;
} complex;

Display* x11setup(Window *win, GC *gc, int width, int height); //Function prototype

int cal_pixel(complex c) {
  int count, max;
  complex z;
  float temp, lengthsq;
  max = 256; //Maximum number of iterations to do
  z.real = 0; z.imaginary = 0; //Initialise value of complex struct
  count = 0; //Current number of iterations
  do {
	temp = z.real * z.real - z.imaginary * z.imaginary + c.real;
    z.imaginary = 2 * z.real * z.imaginary + c.imaginary;
    z.real = temp;
    lengthsq = z.real * z.real + z.imaginary * z.imaginary;
    count++;
  } while ((lengthsq < 4.0) && (count < max));
  return count;
}

int main(int argc, char *argv[])
{
	int rank, worldSize, i, j, x, y;
	double time;
	unsigned int width = X_RESN, height = Y_RESN; //Window size
	Window win; //Initialization for a window
	GC gc; //Graphics context
	Display *display = NULL;
	
	MPI_Init(&argc, &argv); 
	MPI_Comm_size(MPI_COMM_WORLD, &worldSize); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
	if(rank==0) //Master node operations
	{
		int mandelbrot[Y_RESN][X_RESN] = {0}; //2D array to store the mandelbrot pixel values into
		int imageLine, running = 1;
		MPI_Status stat;
		
		display = x11setup(&win, &gc, width, height);
		time = MPI_Wtime(); //Get the start time
		
		for(i=0;i<Y_RESN;++i) { //Recv for the number of times there are lines in the Y resolution
			MPI_Recv(&imageLine, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &stat); //Receive which line will be added
			MPI_Recv(&(mandelbrot[imageLine]), X_RESN, MPI_INT, stat.MPI_SOURCE, 2, MPI_COMM_WORLD, &stat); //Receive mandelbrot line
		}
		
		time = MPI_Wtime() - time; //Get time taken to calculate the mandelbrot
		printf("Calculation time took %f seconds\n", time); //Print elapsed time
		
		XClearWindow(display, win); //Clear window and draw the Mandelbrot
		for(i=0;i<X_RESN;++i) {
			for(j=0;j<Y_RESN;++j) {
				if(mandelbrot[i][j]==256)
					XDrawPoint(display, win, gc, j, i); //Draw point at i,j in white
			//Use XDrawPoint(display, win, gc, x, y) to draw a single pixel at (x,y)
			XFlush(display);
			}
		}

		while(running) { //Wait for user to exit screen with keypress
			if(XPending(display)) {
				XEvent ev;
				XNextEvent(display, &ev);
				switch(ev.type) {
					case KeyPress:
						running = 0;
						break;
				}
			}
		}
		XCloseDisplay(display); //Close the display window
	} //End master node operations
	
	else { //Slave node operations
		complex c;
		int line = rank - 1, nodes = worldSize - 1;
		float realStep = (float)(REAL_MAX - REAL_MIN) / X_RESN, imagStep = (float)(IMAG_MAX - IMAG_MIN) / Y_RESN; //Real and imaginary interpolation values
		int mandelbrotLine[X_RESN] = {0}; //1D array to store line value into
		
		while(line <= Y_RESN) {
			c.real = -2; //Set real and imaginary parts of the complex number to start with
			c.imaginary = 2 - (line * imagStep);
			for(i=0;i<X_RESN;++i) { //Calculate every pixel in the line
				mandelbrotLine[i] = cal_pixel(c);
				c.real += realStep; //Increment the real value of the complex number
			}
			MPI_Send(&line, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
			MPI_Send(&mandelbrotLine, X_RESN, MPI_INT, 0, 2, MPI_COMM_WORLD);
			line += nodes; //Go to the next line
		}
		
	} //End slave node operations
	
	MPI_Finalize();
	return 0;
}


Display * x11setup(Window *win, GC *gc, int width, int height)
{
	
	/* --------------------------- X11 graphics setup ------------------------------ */
	Display 		*display;
	unsigned int 	win_x,win_y, /* window position */
					border_width, /* border width in pixels */
					display_width, display_height, /* size of screen */
					screen; /* which screen */
	
	char 			window_name[] = "Mandelbrot", *display_name = NULL;
	unsigned long 	valuemask = 0;
	XGCValues 		values;
	
	XSizeHints 		size_hints;
	
	//Pixmap 		bitmap;
	//XPoint 		points[800];
	FILE 			*fopen ();//, *fp;
	//char 			str[100];
	
	XSetWindowAttributes attr[1];
	
	if ( (display = XOpenDisplay (display_name)) == NULL ) { /* connect to Xserver */
		fprintf (stderr, "Cannot connect to X server %s\n",XDisplayName (display_name) );
		exit (-1);
	}
	
	screen = DefaultScreen (display); /* get screen size */
	display_width = DisplayWidth (display, screen);
	display_height = DisplayHeight (display, screen);
	
	win_x = 0; win_y = 0; /* set window position */
	
	border_width = 4; /* create opaque window */
	*win = XCreateSimpleWindow (display, RootWindow (display, screen),
			win_x, win_y, width, height, border_width,
			WhitePixel (display, screen), BlackPixel (display, screen));
			
	size_hints.flags = USPosition|USSize;
	size_hints.x = win_x;
	size_hints.y = win_y;
	size_hints.width = width;
	size_hints.height = height;
	size_hints.min_width = 300;
	size_hints.min_height = 300;
	
	XSetNormalHints (display, *win, &size_hints);
	XStoreName(display, *win, window_name);
	
	*gc = XCreateGC (display, *win, valuemask, &values); /* create graphics context */
	
	XSetBackground (display, *gc, BlackPixel (display, screen));
	XSetForeground (display, *gc, WhitePixel (display, screen));
	XSetLineAttributes (display, *gc, 1, LineSolid, CapRound, JoinRound);
	
	attr[0].backing_store = Always;
	attr[0].backing_planes = 1;
	attr[0].backing_pixel = BlackPixel(display, screen);
	
	XChangeWindowAttributes(display, *win, CWBackingStore | CWBackingPlanes | CWBackingPixel, attr);
	
	XSelectInput(display, *win, KeyPressMask);
	
	XMapWindow (display, *win);
	XSync(display, 0);
	
	/* --------------------------- End of X11 graphics setup ------------------------------ */
	return display;
}