
# include <iostream>
# include <cmath>
# include <math.h>
# include <iomanip>
using namespace std;

//////////////////////////////declarations//////////////////////////////

struct vec2;
vec2* vec_ctor(float, float);
void printv(vec2*);
float dot(vec2*, vec2*);
float veclen(vec2*);

struct vec3;
vec3* vec_ctor(float, float, float);
void printv(vec3*);
float dot(vec3*, vec3*);
float veclen(vec3*);

struct mat3x3;
float mat_sum(mat3x3*);
mat3x3* mat_ctor(float[9]);
void mat_write(mat3x3*, float[9]);
void printm(mat3x3*);
void mat_mul(mat3x3*, mat3x3*, mat3x3*);
void mat_scale(mat3x3*, float);

struct cell;
void cell_init(cell*, int, int);

struct grid;
grid* grid_init(int, int);

float sin_beta(int);
float cos_beta(int);
mat3x3* elliptic_mat(vec2*);
mat3x3* diffusion_mat(float, vec2*);

//////////////////////////////2-vector//////////////////////////////

struct vec2 {
	float x;
	float y;

	vec2(float x = 0, float y = 0):x(x), y(y){}
	
	inline vec2& operator = (const vec2& v) {
		x = v.x;
		y = v.y;
		return *this;
	}

	inline vec2 operator + (const vec2& v) const {
		return vec2(v.x + x, v.y + y);
	}
};

vec2* vec_ctor(float x, float y) {
	vec2* v = new vec2;
	v->x = x;
	v->y = y;
	return v;
}

void printv(vec2* v) {
	cout << "( " << v->x << ", " << v->y << " )\n";
}

float dot(vec2* v1, vec2* v2) {
	return (v1->x * v2->x + v1->y * v2->y);
}

float veclen(vec2* v) {
	return sqrt(dot(v, v));
}

//////////////////////////////3-vector//////////////////////////////

struct vec3 {
	float x;
	float y;
	float z;
	
	inline vec3 operator + (vec3 v) {
		return { v.x + x, v.y + y , v.z + z};
	}
};

vec3* vec_ctor(float x, float y, float z) {
	vec3* v = new vec3;
	v->x = x;
	v->y = y;
	v->z = z;
	return v;
}

void printv(vec3* v) {
	cout << "( " << v->x << ", " << v->y << ", " << v->z << " )\n";
}

float dot(vec3* v1, vec3* v2) {
	return (v1->x * v2->x + v1->y * v2->y + v1->z * v2->z);
}

float veclen(vec3* v) {
	return sqrt(dot(v, v));
}

//////////////////////////////3x3 matrix//////////////////////////////

struct mat3x3{
	float nums[3][3] = { 0,0,0,0,0,0,0,0,0 };
};

float mat_sum(mat3x3 *mat) {
	// return sum of all elements
	return (mat->nums[0][0] + mat->nums[0][1] + mat->nums[0][2] + mat->nums[1][0] + mat->nums[1][1] + mat->nums[1][2] + mat->nums[2][0] + mat->nums[2][1] + mat->nums[2][2]);
}

mat3x3* mat_ctor(float nums[9]) {
	mat3x3* mat = new mat3x3;
	for (int i = 0; i < 9; i++) {
		mat->nums[i / 3][i % 3] = nums[i];
	}
	return mat;
}

void mat_write(mat3x3* mat, float nums[9]) {
	for (int i = 0; i < 9; i++) {
		mat->nums[i / 3][i % 3] = nums[i];
	}
}

void printm(mat3x3* mat) {
	cout << fixed;
	cout << "( ";
	for (int i = 0; i < 9; i++) {
		cout << mat->nums[i / 3][i % 3] << " ";
		if ((i % 3 == 2) && (i < 8)) {
			cout << "\n  ";
		}
	}
	cout << " )\n";
}

void mat_mul(mat3x3* out, mat3x3* a, mat3x3* b) {
	float nums[9] = { 0 };
	for (int i = 0; i < 9; i++) {
		nums[i] = a->nums[i / 3][0] * b->nums[0][i % 3] + a->nums[i / 3][1] * b->nums[1][i % 3] + a->nums[i / 3][2] * b->nums[2][i % 3];
	}
	mat_write(out, nums);
}

void mat_scale(mat3x3* mat, float num) {
	for (int i = 0; i < 9; i++) {
		mat->nums[i / 3][i % 3] = num * mat->nums[i / 3][i % 3];
	}
}

//////////////////////////////cell//////////////////////////////

struct cell {
	float temp;
	float ntemp;
	float pressure;
	float npressure;
	vec2* wind;
	float diff_rate;
	mat3x3* diff_mat;
};

void cell_init(cell* c, int x, int y) {
	//TODO generate temp+pressure with perlin
	c = new cell;
	c->temp = 0.0;
	c->ntemp = 0.0;
	c->pressure = 0.0;
	c->npressure = 0.0;
	c->wind = vec_ctor(0.0, 0.0);
	c->diff_rate = 0.0;
	float nums[9] = { 0 };
	c->diff_mat = mat_ctor(nums);
}

//////////////////////////////grid//////////////////////////////

struct grid {
	cell** cells;
	int xmax;
	int ymax;
};

grid* grid_init(int x, int y) {
	grid* g = (grid*)malloc(sizeof(grid));
	if (g == NULL) {
		throw runtime_error("Grid struct allocation fail");
	}
	else {
		g->xmax = x;
		g->ymax = y;
		g->cells = (cell**)malloc(sizeof(cell*) * y);
		if (g->cells == NULL) {
			throw runtime_error("Cell grid allocation fail");
		}
		else {
			for (int i = 0; i < y; i++) {
				g->cells[i] = (cell*)malloc(sizeof(cell) * x);
				if (g->cells[i] == NULL) {
					throw runtime_error("Cell row allocation fail");
				}
				else {
					for (int j = 0; j < x; j++) {
						cell_init(&(g->cells[j][i]), j, i);
					}
				}

			}
		}
		
		return g;
	}
	
}

//////////////////////////////model functions//////////////////////////////

float sin_beta(int dir) {
	float arr[8] = { 0, sqrt(2) / 2, 1, sqrt(2) / 2 , 0, -sqrt(2) / 2 , -1, -sqrt(2) / 2 };
	return arr[dir];
}

float cos_beta(int dir) {
	float arr[8] = { 1, sqrt(2) / 2 , 0, -sqrt(2) / 2 , -1, -sqrt(2) / 2 , 0, sqrt(2) / 2 };
	return arr[dir];
}

mat3x3* elliptic_mat(vec2* wind) {
	const float efd = 1.0;
	// extra-focal distance = semi-major axis - focal distance
	float nums[9] = { 0 };
	int dir_mat_map[8] = { 5, 2, 1, 0, 3, 6, 7, 8 };
	// map from 8 cardinal directions to 3x3 matrix positions
	// 4 is center, nothing maps there
	float wlen = veclen(wind);
	for (int i = 0; i < 8; i++) {
		nums[dir_mat_map[i]] = 2 * efd * (wlen + efd) / (2 * efd + wlen - wind->x * cos_beta(i) - wind->y * sin_beta(i));
	}
	return mat_ctor(nums);
};

mat3x3* diffusion_mat(float diff_rate, vec2* wind) {
	// square normalization matrix
	float snma[9] = { 1,0,0,0,sqrt(2),0,0,0,1 };
	mat3x3* snm = mat_ctor(snma);
	mat3x3* out = elliptic_mat(wind);
	mat_mul(out, snm, out);
	mat_mul(out, out, snm);
	//rescaling
	mat_scale(out, diff_rate / mat_sum(out));
	//adding retention
	out->nums[1][1] = 1 - diff_rate;
	return out;
}

void main() {

	//////////////////////////////TESTING//////////////////////////////
	/*vec2* w = new vec2;
	w->x = 10;
	w->y = 15;*/
	vec2* w = vec_ctor(10, 15);
	printv(w);
	vec2* v = new vec2;
	*v = *w + *w;
	printv(v);
	printm(elliptic_mat(w));
	//float a[9] = { 1,2,3,4,5,6,7,8,9 };
	//float b[9] = { 2,3,5,7,11,13,17,19,23 };
	//mat3x3* test = mat_mul(mat_ctor(a), mat_ctor(b));
	//mat_print(test);
	printm(diffusion_mat(0.5, w));
	////////////////////////////////////////////////////////////
	// 
	//////////////////////////////ACTUAL STUFF//////////////////////////////
	grid* g = grid_init(10, 10);
	//cout << g << "\n";
	//cout << g->xmax << "\n";
	////////////////////////////////////////////////////////////
};