//============================================================================
// Name        :  set_cover.cpp
// Author      : Yamile Patino Vargas
// email			:	ypatino000@citymail.cuny.edu
// Version     : homework 4
// Description : Set Cover Algorithm - Guards
//============================================================================

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <stack>
#include <queue>
#include <stdlib.h>
using namespace std;

FILE *fp;

typedef struct st_orthogonal_polygon{
	XPoint *vertices;
	int  num_vertices;
} orthogonal_polygon;

typedef struct st_covered_sample{
	XPoint *sampleset;
	int  num_samples;
	XPoint candidate;
} covered_sample_set;

orthogonal_polygon *obstacles ;
XPoint *uncovered_samples = NULL;
int num_obstacles = 0;
int uncovered_samples_size = 0;

/*WINDOW DISPLAY: starts*/
Display *display_ptr;
Screen *screen_ptr;
int screen_num;
char *display_name = NULL;
unsigned int display_width, display_height;
Window win;
int border_width, win_x, win_y;
unsigned int win_width, win_height;
XWMHints *wm_hints;
XClassHint *class_hints;
XSizeHints *size_hints;
XTextProperty win_name, icon_name;
char *win_name_string = "Orthogonal Polygons and Guards";
char *icon_name_string;

XEvent report;

GC gc, gc_gray, gc_tred, gc_tyellow, gc_tblue, gc_black;
unsigned long valuemask = 0;
XGCValues  gc_gray_values, gc_black_values, gc_values, gc_tyellow_values, gc_tred_values, gc_tblue_values;
Colormap color_map;
XColor tmp_color1, tmp_color2;

void setColors(Display* display_ptr, Window	win){
	gc = XCreateGC( display_ptr, win, valuemask, &gc_values);
	XSetForeground( display_ptr, gc, BlackPixel( display_ptr, screen_num ) );
	XSetLineAttributes( display_ptr, gc, 4, LineSolid, CapRound, JoinRound);

	/* and three other graphics contexts, to draw in yellow and red and grey white=WhitePixel(dis, screen);*/
	gc_black = XCreateGC( display_ptr, win, valuemask, &gc_black_values);
	XSetLineAttributes(display_ptr, gc_black, 1, LineSolid,CapRound, JoinRound);
	if( XAllocNamedColor( display_ptr, color_map, "black",
			&tmp_color1, &tmp_color2 ) == 0 )
	{printf("failed to get color black\n"); exit(-1);}
	else
		XSetForeground( display_ptr, gc_black, BlackPixel(display_ptr, screen_num ) );

	gc_tyellow = XCreateGC( display_ptr, win, valuemask, &gc_tyellow_values);
	XSetLineAttributes(display_ptr, gc_tyellow, 1, LineSolid,CapRound, JoinRound);
	if( XAllocNamedColor( display_ptr, color_map, "yellow", &tmp_color1, &tmp_color2 ) == 0 )
	{printf("failed to get color yellow\n"); exit(-1);}
	else{
		XSetForeground( display_ptr, gc_tyellow, tmp_color1.pixel );
	}

	/* other graphics contexts red*/
	gc_tred = XCreateGC( display_ptr, win, valuemask, &gc_tred_values);
	XSetLineAttributes( display_ptr, gc_tred, 3, LineSolid, CapRound, JoinRound);
	if( XAllocNamedColor( display_ptr, color_map, "red", &tmp_color1, &tmp_color2 ) == 0 )
	{printf("failed to get color red\n"); exit(-1);}
	else
	{XSetForeground( display_ptr, gc_tred, tmp_color1.pixel );}
	/* other graphics contexts red*/

	gc_gray = XCreateGC( display_ptr, win, valuemask, &gc_gray_values);
	XSetLineAttributes( display_ptr, gc_gray, 2, LineSolid, CapRound, JoinRound);
	if( XAllocNamedColor( display_ptr, color_map, "light grey", &tmp_color1, &tmp_color2 ) == 0 )
	{printf("failed to get color gray\n"); exit(-1);}
	else
		XSetForeground( display_ptr, gc_gray, tmp_color1.pixel );

	/*other graphics contexts grey*/
	gc_tblue = XCreateGC( display_ptr, win, valuemask, &gc_tblue_values);
	XSetLineAttributes( display_ptr, gc_tblue, 1, LineSolid, CapRound, JoinRound);
	if( XAllocNamedColor( display_ptr, color_map, "blue", &tmp_color1, &tmp_color2 ) == 0 )
	{printf("failed to get color blue\n"); exit(-1);}
	else{
		XSetForeground( display_ptr, gc_tblue, tmp_color1.pixel );
	}
}

/*Distance between two points*/
double distance_p(XPoint p1, XPoint p2){
	double l = 0;
	if(!(p1.x==p2.x and p1.y == p2.y))
		l = (sqrt( pow( p2.x - p1.x, 2 ) + pow(p2.y - p1.y, 2 )));
	return l;
}

XPoint  getcoordinates(char* coordinates){
	char* x=0;
	char* y=0;
	char* temp;
	XPoint p;

	temp = strtok(coordinates, ",");

	while (temp)
	{//get cord_x, cord_y per every strtok
		char* open_par = strchr(temp, '(');
		if(open_par){
			*open_par = 0;
			open_par++;
			x = new char[(strlen(open_par) + 1)];
			strcpy(x, open_par);
			p.x = (short)strtol(x, NULL, 10);
		}else{
			char* close_par = strchr(temp, ')');
			*close_par=0;
			y = new char[(strlen(close_par) + 1)];
			strcpy(y, temp);
			p.y = (short)strtol(y, NULL, 10);
		}
		temp = strtok(0, ",");
	}
	free(x);  free(y);  free(temp);
	return p;
}

orthogonal_polygon getvertex(char* line){
	orthogonal_polygon polygon;
	polygon.num_vertices =0;
	char* s;
	int num_vertex=0;
	int lenght_s = (strlen(line) /4);
	if(lenght_s !=0){
		char **points = new char*[ (strlen(line) /4)];
		s = strtok(line, " ");
		while (s){
			if(strlen(s) >= 5){
			points[num_vertex] = new char[(strlen(s) + 1)];
			strcpy(points[num_vertex], s);
			num_vertex++;
			}
			s = strtok(0, " ");
		}
		XPoint *vertex = new XPoint[num_vertex];
		//get cord_x, cord_y per every strtok
		int i = 0; bool valid_polygon = true;
		vertex[num_vertex-1] = getcoordinates( points[num_vertex-1]);
		XPoint temp_vertex = vertex[num_vertex-1] ;
		for (i=0; i< num_vertex-1; i++){
			vertex[i] = getcoordinates( points[i]);
			if(vertex[i].x != temp_vertex.x and vertex[i].y != temp_vertex.y){
				valid_polygon = false;
				i=num_vertex;
			}else
				temp_vertex = vertex[i];
		}
		if(valid_polygon){
			polygon.num_vertices = num_vertex;
			polygon.vertices = vertex;
		}
	}
	return polygon;
}

void read_file(FILE *t){
	stack <orthogonal_polygon>  the_obstacles;
	char *temp_line  = new char[1000];
	orthogonal_polygon temp_pol;
	int i =0;
	while(!feof(t)){
		if ( fgets (temp_line ,1000 , t) != NULL ){
			temp_pol = getvertex(temp_line);
			if(temp_pol.num_vertices != 0){
				the_obstacles.push( temp_pol);
				i++;
			}
			else
				cout << "Invalid polygon or Empty line at # " << i << endl;
		}
	}
	fclose(t);
	num_obstacles = the_obstacles.size();
	obstacles = new orthogonal_polygon[num_obstacles];
	for(i=0; !the_obstacles.empty(); i++){
		obstacles[i] = the_obstacles.top();
		the_obstacles.pop();
		XFillPolygon(display_ptr, win, gc_gray,obstacles[i].vertices,obstacles[i].num_vertices , Nonconvex, CoordModeOrigin);
		XDrawLines(display_ptr, win, gc_black,obstacles[i].vertices,obstacles[i].num_vertices ,  CoordModeOrigin);
		XDrawLine(display_ptr, win, gc_black,obstacles[i].vertices[0].x, obstacles[i].vertices[0].y,
				obstacles[i].vertices[obstacles[i].num_vertices-1].x,  obstacles[i].vertices[obstacles[i].num_vertices-1].y);
	}
	free(temp_line);
}

/*
 * The merge_sort algorithm sorts the array of subsets by the increasing number of ones (bites)
 * */
covered_sample_set *merge_sort(covered_sample_set *sample_sets, int n){
	int i = 0;
	covered_sample_set * left_list = NULL;
	covered_sample_set * right_list = NULL;
	if(n == 1)
		return sample_sets;
	left_list = new covered_sample_set[n/2];
	right_list = new covered_sample_set[n-(n/2)];

	for (i = 0; i< n/2;i++){
		left_list[i]=sample_sets[i];
		right_list[i]=sample_sets[i+(n/2)];
	}
	if(n%2!=0)
		right_list[n/2]=sample_sets[n-1];

	/*Recursion*/
	right_list = merge_sort(right_list, (n - n/2));
	left_list = merge_sort(left_list, n/2);

	covered_sample_set *final_list = new covered_sample_set[n];

	i = 0; int j = 0; int fin = 0;
	while((i < n/2) and (j < (n - n/2))){
		if(left_list[i].num_samples  >= right_list[j].num_samples){
			final_list[fin] = left_list[i];
			fin++; i++;
		}else{
			final_list[fin] = right_list[j];
			fin++; j++;
		}
	}
	if(i < n/2){
		while(i < n/2){
			final_list[fin] = left_list[i];
			fin++; i++;
		}
	}else if(j < (n - n/2)){
		while(j < (n - n/2)){
			final_list[fin] = right_list[j];
			fin++; j++;
		}
	}
	left_list = 0;
	right_list = 0;
	return final_list;
}


void set_uncovered_size(int s){
	uncovered_samples_size = s;
}
void set_uncovered_samples(XPoint *S){
	uncovered_samples = S;
}

XPoint *get_uncovered_samples(){
	return uncovered_samples;
}
int get_uncovered_size(){
	return uncovered_samples_size;
}
void draw_guards_covered_samples(stack <covered_sample_set> guards){
	int i =0;
	int num_samples = 0;
	while(!guards.empty()){
	num_samples = guards.top().num_samples;
		for(i=0; i< num_samples; i++){
			XDrawLine(display_ptr, win, gc_tyellow, guards.top().candidate.x, guards.top().candidate.y, guards.top().sampleset[i].x, guards.top().sampleset[i].y);
			XDrawRectangle(display_ptr, win, gc_tblue, guards.top().sampleset[i].x, guards.top().sampleset[i].y, 1,1);
		}
		XDrawRectangle(display_ptr, win, gc_tred, guards.top().candidate.x, guards.top().candidate.y, 1,1);
		guards.pop();
	}
}
void draw_uncovered(XPoint *uncovered, int size){
	int i=0;
	for(i=0; i < size; i++){
		XDrawRectangle(display_ptr, win, gc_black, uncovered[i].x, uncovered[i].y, 1,1);
	}
}

/*intersect: Evaluate if the points intersect.
 * Existent edge from p to q.
 * Temp edge from r to s
 * returns true if the edge pq separates rs.
 * and false otherwise
 */
bool intersect(XPoint s, XPoint r, XPoint p, XPoint q){
	double   orientations_pqr, orientations_pqs, orientations_rsp, orientations_rsq, result_pq, result_rs;
	int px,py,qx,qy,sx,sy,rx,ry;

	px=p.x; py=p.y; qx=q.x; qy=q.y; sx=s.x; sy=s.y; rx=r.x; ry=r.y;
	orientations_pqr = ((px * qy) + (py * rx) + (qx * ry)) - ((px * ry) + (qx * py) + (rx * qy) );
	orientations_pqs = ((px * qy) + (py * sx) + (qx * sy)) - ((px * sy) + (qx * py) + (sx * qy) );
	result_pq = orientations_pqr * orientations_pqs;
	orientations_rsp = ((rx * sy) + (ry * px) + (sx * py)) - ((rx * py) + (sx * ry) + (px * sy) );
	orientations_rsq = ((rx * sy) + (ry * qx) + (sx * qy)) - ((rx * qy) + (sx * ry) + (qx * sy) );
	result_rs = orientations_rsp * orientations_rsq;

	return( ((result_pq < 0) and (result_rs < 0)) or (orientations_rsp == 0) or (orientations_rsq == 0));
}

int evaluate_point(XPoint p, orthogonal_polygon obstacle){
	int count_intersections=0; int i =0;
	XPoint temp_vertex = obstacle.vertices[obstacle.num_vertices-1] ;
	for(i=0; i< obstacle.num_vertices;i++){
		if((p.y > min(temp_vertex.y, obstacle.vertices[i].y)) and (p.y <= max(temp_vertex.y, obstacle.vertices[i].y)) and (p.x <= temp_vertex.x) ){
			if(p.x == temp_vertex.x){
				i= obstacle.num_vertices;
				count_intersections = 2;
			}else{
				count_intersections++;
			}
		}
		temp_vertex = obstacle.vertices[i];
	}
	return count_intersections;
}

int evaluate_candidate(XPoint p, orthogonal_polygon obstacle){
	int count_intersections=0; int i =0;
	XPoint temp_vertex = obstacle.vertices[obstacle.num_vertices-1] ;
	for(i=0; i< obstacle.num_vertices;i++){
		if((p.y > min(temp_vertex.y, obstacle.vertices[i].y)) and (p.y <= max(temp_vertex.y, obstacle.vertices[i].y)) and (p.x <= temp_vertex.x) ){
			if(p.x == temp_vertex.x  ){
				i= obstacle.num_vertices;
				count_intersections = 1;
			}else
				count_intersections++;
		}else if(p.y == temp_vertex.y and p.y == temp_vertex.y){
			i= obstacle.num_vertices;
			count_intersections = 1;
		}
		temp_vertex = obstacle.vertices[i];
	}
	return count_intersections;
}
/*Per evry obstacle evaluate if the point is inside or out, when you find that there is inside one the function terminates */
bool isvalidapoint(XPoint p, orthogonal_polygon *obstacle_s){
	int i=0;
	bool valid = true;
	for (i=0; i<num_obstacles; i++){
		if(evaluate_point(p,obstacle_s[i])%2 != 0){
			valid = false;
			i= num_obstacles;
		}
	}
	return valid;
}
bool isvalidcandidate(XPoint p, orthogonal_polygon *obstacle_s){
	//Per evry obstacle evaluate if the point is inside or out, when you find that there is inside one then program terminates
	int i=0;
	bool valid = true;
	for (i=0; i<num_obstacles; i++){
		if(evaluate_candidate(p,obstacle_s[i])%2 != 0){
			valid = false;
			i= num_obstacles;
		}
	}
	return valid;
}
/*  samplepoints() returns a dense set of valid points (point outside a obstacle), every obstacle must have
 *   points around its sides
 * */
XPoint *samplepoints(orthogonal_polygon *obstacle_s, int size){
	XPoint sp; XPoint *sample_points = new XPoint[size];
	int x,y,valid_points,i;
	valid_points = i = x= y =0;
	while(valid_points< size){
		sp.x = (rand() % (rand() % 42389 ) ) % 500;
		sp.y = (rand() % (rand() % 18756 + 879) ) % 500;

		if (isvalidapoint(sp, obstacle_s)){
			sample_points[valid_points] = sp;
			valid_points++;
		}
	}
	return sample_points;
}

XPoint *candidatepoints(orthogonal_polygon *obstacle_s, int size){
	XPoint sp; XPoint *candidate_points = new XPoint[size];
	//generate random points on the plane [0, 499] and evaluate if they are valid.
	int x,y,valid_points,i;
	valid_points = i = x = y =0;
	while(valid_points< size){
		sp.x = (rand() % 30000 ) % 500;
		sp.y = (rand() % 25000 +500) % 500;

		if (isvalidcandidate(sp, obstacle_s)){
			candidate_points[valid_points] = sp;
			valid_points++;
		}
	}
	return candidate_points;
}
/*
 * Iintersects_obstacles returns true when the edge c to s intersect obstacles.
 * */
bool intersects_obstacles(XPoint c, XPoint s, orthogonal_polygon *obs){
	int i=0;int j=0; int n = 0 ;
	XPoint temp;
	bool exist_intersection = false;
	for(i=0; i < num_obstacles; i++){
		temp = obs[i].vertices[obs[i].num_vertices-1];
		n = obs[i].num_vertices;
		for(j=0; j < n; j++){
			if(intersect( c, s, obs[i].vertices[j], temp)){
				exist_intersection = true;
				j = obs[i].num_vertices;
				i =  num_obstacles;
			}else
				temp = obs[i].vertices[j];
		}
	}
	return exist_intersection;
}

/*Returns all sample points visible from c
 * */
covered_sample_set find_covered_samples(XPoint c, XPoint *s, int nums){
	covered_sample_set temp_css;
	stack <XPoint> coveredsamples;
	int i = 0;
	//if c s[i] does not intersect any edge of the obstacles put it in stack
	for(i=0; i<nums; i++){
		if(!intersects_obstacles(c,s[i], obstacles)){
			coveredsamples.push(s[i]);
		}
	}
	temp_css.num_samples = coveredsamples.size();
	temp_css.sampleset = new XPoint[temp_css.num_samples ];

	for(i=0; i< temp_css.num_samples ; i++){
		temp_css.sampleset[i] = coveredsamples.top();
		coveredsamples.pop();
	}
	temp_css.candidate = c;
	return temp_css;
}

bool belongsto(XPoint s, covered_sample_set c){
	bool belongs = false;
	int i=0;
	for(i=0; i< c.num_samples; i++){
		if(s.x == c.sampleset[i].x and s.y == c.sampleset[i].y ){
			belongs = true;
			i = c.num_samples;
		}
	}
	return belongs;
}

XPoint *rest_samples(covered_sample_set best_candidate_set,  XPoint *S, int nums){
	XPoint *new_s;
	stack <XPoint> temp_s;
	int i=0;
	for(i =0; i < nums; i++){
		if(!belongsto(S[i], best_candidate_set)){
			temp_s.push(S[i]);
		}
	}
	int j = temp_s.size();
	set_uncovered_size(temp_s.size());
	new_s = new XPoint[j];
	for(i=0; i< j; i++){
		new_s[i] = temp_s.top();
		temp_s.pop();
	}
	return new_s;
}

XPoint *rest_candidates(XPoint c,  XPoint *candidates, int numc){
	XPoint *rest_cand = new XPoint[numc - 1];
	int i=0; int j= 0;
	for(i=0; i< numc ; i++){
		if(c.x == candidates[i].x and c.y == candidates[i].y){
			j=i+1;
			i = numc;
		}else{
			rest_cand[i]= candidates[i];
		}
	}
	while(j <numc ){
		rest_cand[j-1]= candidates[j];
		j++;
	}
	return rest_cand;
}

stack<covered_sample_set> greedy_set_cover_aproximation(XPoint *C, int numc,  XPoint *S, int nums, int orignal_nums){
	stack <covered_sample_set> candidates;
	XPoint *S_s = S;
	XPoint *C_c = C;
	double cond = orignal_nums*0.05 < 1? 1: orignal_nums*0.05;
	if(nums <= cond or numc == 0){
		set_uncovered_samples(S_s);
		set_uncovered_size(nums);
		return candidates;
	}
	covered_sample_set *css = new covered_sample_set[numc];
	int i,j;
	i = j = 0;
	for(i = 0; i<numc; i++){
		css[i] = find_covered_samples(C_c[i], S_s, nums);
	}
	css = merge_sort(css, numc);
	if(css[0].num_samples != 0){
		S_s = rest_samples(css[0], S, nums);
		C_c = rest_candidates(css[0].candidate, C , numc);
		candidates = greedy_set_cover_aproximation(C_c, numc - 1,  S_s, get_uncovered_size(), orignal_nums) ;
		candidates.push(css[0]);
	}else{
		set_uncovered_samples(S_s);
		set_uncovered_size(nums);
	}
	css = 0;
	return candidates;
}

void findguards(){
	int num_samples = 1000;
	int num_candidates = 1000;
	XPoint *S = samplepoints(obstacles, num_samples);
	XPoint  *C = candidatepoints(obstacles, num_candidates);
	stack <covered_sample_set> guard_samples = greedy_set_cover_aproximation(C, num_candidates,  S, num_samples, num_samples) ;

	double uncovered_perct =   (double)(get_uncovered_size() * 100)/num_samples;
	double covered_perct =  (double) ((num_samples-get_uncovered_size() )* 100)/num_samples;
	cout << "Sample set: " << endl;
	cout << guard_samples.size() << " Guards cover  " << covered_perct  << "% out of " << num_samples << " sample points."<< endl;
	cout << uncovered_perct  << "% of the samples were not covered. "<< endl;

	while(uncovered_perct > 5 and num_samples < 4000){
		cout << "New sample set: " << endl;
		set_uncovered_samples(0);
		set_uncovered_size(num_samples);
		num_samples = num_samples * 1.2;
		num_candidates = num_candidates * 1.2;
		S=0; C=0;
		S = samplepoints(obstacles, num_samples);
		C = candidatepoints(obstacles, num_candidates);
		guard_samples = greedy_set_cover_aproximation(C, num_candidates,  S, num_samples, num_samples) ;
		uncovered_perct =    (double)(get_uncovered_size() * 100)/num_samples;
		covered_perct =  (double) ((num_samples-get_uncovered_size() )* 100)/num_samples;

		cout << guard_samples.size() << " Guards cover  " << covered_perct << "% out of " << num_samples << " sample points."<< endl;
		cout << uncovered_perct  << "% of the samples were not covered. "<< endl;cout << endl;
	}
	draw_guards_covered_samples(guard_samples);
	draw_uncovered( get_uncovered_samples(), get_uncovered_size());

	C = 0;
	S = 0;
}

int main(int argc, char *argv[])
{
	/* opening display: basic connection to X Server */
	if( (display_ptr = XOpenDisplay(display_name)) == NULL ){printf("Could not open display. \n"); exit(-1);}
	screen_num = DefaultScreen( display_ptr );
	screen_ptr = DefaultScreenOfDisplay( display_ptr );
	color_map  = XDefaultColormap( display_ptr, screen_num );
	display_width  = DisplayWidth( display_ptr, screen_num );
	display_height = DisplayHeight( display_ptr, screen_num );
	/* creating the window */
	border_width = 10;
	win_x = 0; win_y = 0; win_width = 520; win_height = 510;
	/*rectangular window*/

	win= XCreateSimpleWindow( display_ptr, RootWindow( display_ptr, screen_num),
			win_x, win_y, win_width, win_height, border_width,
			BlackPixel(display_ptr, screen_num),
			WhitePixel(display_ptr, screen_num) );

	/* now try to put it on screen, this needs cooperation of window manager */
	size_hints = XAllocSizeHints();
	wm_hints = XAllocWMHints();
	class_hints = XAllocClassHint();
	if( size_hints == NULL || wm_hints == NULL || class_hints == NULL )
	{ printf("Error allocating memory for hints. \n"); exit(-1);}

	size_hints -> flags = PPosition | PSize | PMinSize  ;
	size_hints -> min_width = 60;
	size_hints -> min_height = 60;

	XStringListToTextProperty( &win_name_string,1,&win_name);
	XStringListToTextProperty( &icon_name_string,1,&icon_name);

	wm_hints -> flags = StateHint | InputHint ;
	wm_hints -> initial_state = NormalState;
	wm_hints -> input = False;

	XSetWMProperties( display_ptr, win, &win_name, &icon_name, argv, argc, size_hints, wm_hints, class_hints );

	/* what events do we want to receive */
	XSelectInput( display_ptr, win, ExposureMask | StructureNotifyMask | ButtonPressMask );

	/* finally: put window on screen */
	XMapWindow( display_ptr, win );
	XFlush(display_ptr);

	/* create graphics context, so that we may draw in this window */
	setColors(display_ptr,win);
	if(argc > 1){
		fp = fopen(argv[1], "r");
		read_file(fp);
	}else
		printf("File was not given.\n");
	findguards();
	while(1)
	{
		XNextEvent( display_ptr, &report );
		switch( report.type )
		{
		case Expose:
			break;
		case ConfigureNotify:
			/* This event happens when the user changes the size of the window*/
			win_width = report.xconfigure.width;
			win_height = report.xconfigure.height;
			break;
		case ButtonPress:
		{
		if (report.xbutton.button == Button3){
				XFreeGC(display_ptr, gc);
				XFreeGC(display_ptr, gc_gray);
				XFreeGC(display_ptr, gc_tred);
				XDestroyWindow(display_ptr,win);
				XCloseDisplay(display_ptr);
			}
		}
		break;
		default:
			/* this is a catch-all for other events; it does not do anything.
	              One could look at the report type to see what the event was */
			break;
		}
	}
	exit(0);
}
