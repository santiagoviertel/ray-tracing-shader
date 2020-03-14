// Author:	Santiago Viertel
// Title:	Ray Tracing Shader

precision mediump float;

uniform vec2 u_resolution;
uniform float u_time;

#define HEAPSIZE		7				// Number of nodes of the heap
#define EPSILON			0.000001		// Small constant for avoid overflow
#define AIRRFI			1.00027316		// Refractive index of the air
const vec3 BACKGROUND	= vec3(1.0);	// Background color

struct camera {
	vec3 pos;	// Position
	vec3 upv;	// Up vector
	vec3 lav;	// Look at vector
	vec3 siv;	// Side vector
};

struct lightComponents {
	vec3 dif;	// Diffuse component
	vec3 spe;	// Specular component
};

struct pointLight {
	vec3 pos;				// Position
	lightComponents com;	// Light components
};

struct material {
	vec3 amb;	// Ambient component
	vec3 dif;	// Diffuse component
	vec3 spe;	// Specular component
	float shi;	// Shininess component
	float lap;	// Light absorption percentage
	float rfi;	// Refractive index
};

struct referenceSystem {
	vec3 o;	// Origin of the reference system
	mat3 v;	// The three vectors of the reference system
};

struct rhombusBox {
	int id;				// Object id number
	float halfSide;		// Half of the side of the rhombus/box
	referenceSystem rs;	// The reference system of the rhombus/box. Z is the normal vector
	material mat;		// Material
};

struct sphere {
	int id;			// Object id number
	float rad;		// Radius
	vec3 pos;		// Position
	material mat;	// Material
};

struct ray {
	vec3 o;	// Ray origin
	vec3 d;	// Ray direction
};

struct rhombusSphereDistanceReturn {
	bool i;		// Intersects or do not
	float d;	// Distance from the origin of the ray
};

struct boxDistanceReturn {
	bool i;		// Intersects or do not
	float d;	// Distance from the origin of the ray
	vec3 n;		// Normal vector of the surface
};

struct detectObjectReturn {
	bool inter;		// Intersects or do not
	int idObj;		// Id number of the intersecting object
	vec3 pointObj;	// Point in the object surface
	vec3 normObj;	// Normal vector of the point in the object surface
};

struct heapNode {
	bool exists;	// The node existence
	int idObj;		// Id number of the object from where the ray leave
	ray r;			// The ray parameters
	float mult;		// The multiplier of the color
	vec3 color;		// The color of the node
};

// Simulated objects
camera cam;
float ambientLight;
pointLight lig[2];
rhombusBox floor;
sphere sph;
rhombusBox cub[2];

// Camera rotation
void rotateCamera() {
	vec3 c = 0.5*(sph.pos + cub[0].rs.o + cub[1].rs.o);
	float t = c[1]/cam.lav[1];
	vec3 p1 = cam.pos + t*cam.lav;
	t = c[1]/cam.upv[1];
	vec3 p2 = cam.pos + t*cam.upv;
	t = length(vec2(cam.pos[0], cam.pos[1]) - vec2(p1[0], p1[1]));
	float ang = 0.25*u_time;
	cam.pos[0] = p1[0] + t*sin(ang);
	cam.pos[1] = p1[1] + t*cos(ang);
	cam.lav = normalize(p1 - cam.pos);
	cam.upv = normalize(p2 - cam.pos);
	cam.siv = cross(cam.lav, cam.upv);
}

// Simulated objects initialization
void init() {
	// Camera initialization
	cam = camera(
		vec3( 1.96521,-0.500505, 0.327547),		// Position
		vec3(-0.2834438, 0.1854804, 0.9408808),	// Up vector
		vec3(-0.7872955, 0.5151917,-0.3387379),	// Look at vector
		vec3( 0.5475632, 0.8367643, 0.0)		// Side vector
	);
	// Lights initialization
	ambientLight = 1.0;
	lig[0] = pointLight(
		vec3( 0.027962,-0.218781, 3.28725),	// Position
		lightComponents(
			vec3(1.0),						// Diffuse component
			vec3(1.0)						// Specular component
		)
	);
	lig[1] = pointLight(
		vec3(-3.2987, 2.66388, 3.49542),	// Position
		lightComponents(
			vec3(1.0),						// Diffuse component
			vec3(1.0)						// Specular component
		)
	);
	// Materials initialization
	material ruby = material(
		vec3( 0.1745, 0.01175, 0.01175),		// Ambient component
		vec3( 0.61424, 0.04136, 0.04136),		// Diffuse component
		vec3( 0.727811, 0.626959, 0.626959),	// Specular component
		76.8,									// Shininess component
		0.55,									// Light absorption percentage
		1.76575									// Refractive index
	);
	material emerald = material(
		vec3( 0.0215, 0.1745, 0.0215),		// Ambient component
		vec3( 0.07568, 0.61424, 0.07568),	// Diffuse component
		vec3( 0.633, 0.727811, 0.633),		// Specular component
		76.8,								// Shininess component
		0.55,								// Light absorption percentage
		1.58225								// Refractive index
	);
	material saphire = material(
		vec3( 0.019, 0.042, 0.100),		// Ambient component
		vec3( 0.041, 0.047, 0.510),		// Diffuse component
		vec3( 0.633, 0.633, 0.727811),	// Specular component
		76.8,							// Shininess component
		0.55,							// Light absorption percentage
		1.76575							// Refractive index
	);
	material silver = material(
		vec3( 0.23125, 0.23125, 0.23125),		// Ambient component
		vec3( 0.2775, 0.2775, 0.2775),			// Diffuse component
		vec3( 0.773911, 0.773911, 0.773911),	// Specular component
		89.6,									// Shininess component
		1.0,									// Light absorption percentage
		0.24437									// Refractive index
	);
	// Objects initialization
	floor = rhombusBox(
		1,									// Object id number
		4.0,								// Half of the side of the rhombus
		referenceSystem(
			vec3(-1.8, 1.75751,-1.50001),	// Origin of the reference system
			mat3(
				1.0, 0.0, 0.0,				// X vector
				0.0, 1.0, 0.0,				// Y vector
				0.0, 0.0, 1.0				// Z vector (normal)
			)
		),
		silver							// Material
	);
	sph = sphere(
		2,								// Object id number
		0.590,							// Radius
		vec3(-1.45325, 2.88018,-0.91),	// Position
		ruby							// Material
	);
	cub[0] = rhombusBox(
		3,									// Object id number
		0.528,								// Half of the side of the box
		referenceSystem(
			vec3(-2.7835, 1.54134,-0.972),	// Origin of the reference system
			mat3(
				1.0, 0.0, 0.0,				// X vector
				0.0, 1.0, 0.0,				// Y vector
				0.0, 0.0, 1.0				// Z vector
			)
		),
		emerald								// Material
	);
	cub[1] = rhombusBox(
		4,										// Object id number
		0.435,									// Half of the side of the box
		referenceSystem(
			vec3(-0.893145, 1.24539,-1.065),	// Origin of the reference system
			mat3(
				 0.6678326, 0.7443116, 0.0,		// X vector
				-0.7443116, 0.6678326, 0.0,		// Y vector
				 0.0, 0.0, 1.0					// Z vector
			)
		),
		saphire									// Material
	);
	rotateCamera();
}

// Method of intersection detection between ray and rhombus presented in the paper:
// http://graphics.cs.kuleuven.be/publications/LD04ERQIT/LD04ERQIT_paper.pdf
rhombusSphereDistanceReturn rhombusDistance(ray r, rhombusBox ro, int idObj) {
	if(idObj==ro.id)
		return rhombusSphereDistanceReturn(false, 0.0);
	vec3 v00 = ro.rs.o + ro.halfSide*(ro.rs.v[0] + ro.rs.v[1]);
	vec3 v10 = ro.rs.o + ro.halfSide*(ro.rs.v[0] - ro.rs.v[1]);
	vec3 v01 = ro.rs.o + ro.halfSide*(ro.rs.v[1] - ro.rs.v[0]);
	vec3 e01 = v10 - v00;
	vec3 e03 = v01 - v00;
	vec3 vect1 = cross(r.d, e03);
	float det = dot(e01, vect1);
	if(abs(det) < EPSILON)
		return rhombusSphereDistanceReturn(false, 0.0);
	vec3 vect2 = r.o - v00;
	float alpha = dot(vect2, vect1)/det;
	if(alpha<0.0 || alpha>1.0)
		return rhombusSphereDistanceReturn(false, 0.0);
	vec3 vect3 = cross(vect2, e01);
	float beta = dot(r.d, vect3)/det;
	if(beta<0.0 || beta>1.0)
		return rhombusSphereDistanceReturn(false, 0.0);
	float t = dot(e03, vect3)/det;
	if(t <= 0.0)
		return rhombusSphereDistanceReturn(false, 0.0);
	return rhombusSphereDistanceReturn(true, t);
}

// Method of intersection detection between ray and box
// presented in the book Real-Time Rendering - Third Edition
boxDistanceReturn boxDistance(ray r, rhombusBox bo, int idObj) {
	bool flagMin = true;
	bool flagMax = true;
	float tMin, tMax;
	vec3 nMin, nMax;
	vec3 p = bo.rs.o - r.o;
	for(int i=0;i<3;++i) {
		float e = dot(bo.rs.v[i], p);
		float f = dot(bo.rs.v[i], r.d);
		if(abs(f) > EPSILON) {
			float t1 = (e + bo.halfSide)/f;
			float t2 = (e - bo.halfSide)/f;
			if(t1 < t2) {
				if(flagMin || t1>tMin) {
					tMin = t1;
					nMin = bo.rs.v[i];
					flagMin = false;
				}
				if(flagMax || t2<tMax) {
					tMax = t2;
					nMax = -bo.rs.v[i];
					flagMax = false;
				}
			} else {
				if(flagMin || t2>tMin) {
					tMin = t2;
					nMin = -bo.rs.v[i];
					flagMin = false;
				}
				if(flagMax || t1<tMax) {
					tMax = t1;
					nMax = bo.rs.v[i];
					flagMax = false;
				}
			}
			if(tMin>tMax || tMax<0.0)
				return boxDistanceReturn(false, 0.0, vec3(0.0));
		} else if(-e-bo.halfSide>0.0 || bo.halfSide-e<0.0)
			return boxDistanceReturn(false, 0.0, vec3(0.0));
	}
	if(idObj == bo.id)
		if(abs(tMin) > abs(tMax))
			return boxDistanceReturn(false, 0.0, vec3(0.0));
		else
			return boxDistanceReturn(true, tMax, nMax);
	return boxDistanceReturn(true, tMin, nMin);
}

// Method of intersection detection between ray and sphere
// presented in the book Real-Time Rendering - Third Edition
rhombusSphereDistanceReturn sphereDistance(ray r, sphere sp, int idObj) {
	vec3 l = sp.pos - r.o;
	float s1 = dot(l, r.d);
	float ls = dot(l, l);
	float rs = sp.rad*sp.rad;
	float ms = ls - s1*s1;
	if((ls>=rs && s1<=0.0) || ms>=rs)
		// Origin of the ray after the sphere center AND outside the sphere OR
		// ray above or below the sphere
		return rhombusSphereDistanceReturn(false, 0.0);
	float q = sqrt(rs - ms);
	float t1 = s1 - q;
	float t2 = s1 + q;
	if(idObj == sp.id)
		if(abs(t1) > abs(t2))
			return rhombusSphereDistanceReturn(false, 0.0);
		else
			return rhombusSphereDistanceReturn(true, t2);
	return rhombusSphereDistanceReturn(true, t1);
}

// Method that finds an object that intersects with the ray
// idObj identifies the object from where the ray comes out
detectObjectReturn detectObject(ray r, int idObj) {
	int id = 0;
	float d;
	// Floor
	rhombusSphereDistanceReturn retRhombusSphere = rhombusDistance(r, floor, idObj);
	if(retRhombusSphere.i) {
		id = floor.id;
		d = retRhombusSphere.d;
	}
	// Sphere
	retRhombusSphere = sphereDistance(r, sph, idObj);
	if(retRhombusSphere.i && (id==0 || retRhombusSphere.d<d)) {
		id = sph.id;
		d = retRhombusSphere.d;
	}
	vec3 n;
	// Boxes
	for(int i=0;i<2;++i) {
		boxDistanceReturn retBox = boxDistance(r, cub[i], idObj);
		if(retBox.i && (id==0 || retBox.d<d)) {
			id = cub[i].id;
			n = retBox.n;
			d = retBox.d;
		}
	}
	if(id == 0)
		return detectObjectReturn(false, 0, vec3(0.0), vec3(0.0));
	vec3 p = r.o + d*r.d;
	if(id == floor.id)
		n = floor.rs.v[2];
	else if(id == sph.id)
		n = normalize(p - sph.pos);
	return detectObjectReturn(true, id, p, n);
}

// Light components that the point receive
lightComponents lightReceived(pointLight lig, vec3 ligDir, int idObj, vec3 pointObj) {
	ray r = ray(pointObj, ligDir);
	lightComponents ret = lig.com;
	vec3 pointToLight = lig.pos - pointObj;
	float sDistPointToLight = dot(pointToLight, pointToLight);
	if(idObj != floor.id) {
		rhombusSphereDistanceReturn retRhombus = rhombusDistance(r, floor, idObj);
		if(retRhombus.i && retRhombus.d*retRhombus.d<sDistPointToLight) {
			ret.dif *= (1.0 - floor.mat.lap)*floor.mat.dif;
			ret.spe *= (1.0 - floor.mat.lap)*floor.mat.spe;
		}
	}
	if(idObj != sph.id) {
		rhombusSphereDistanceReturn retSphere = sphereDistance(r, sph, idObj);
		if(retSphere.i && retSphere.d*retSphere.d<sDistPointToLight) {
			ret.dif *= (1.0 - sph.mat.lap)*sph.mat.dif;
			ret.spe *= (1.0 - sph.mat.lap)*sph.mat.spe;
		}
	}
	for(int i=0;i<2;++i)
		if(idObj != cub[i].id) {
			boxDistanceReturn retBox = boxDistance(r, cub[i], idObj);
			if(retBox.i && retBox.d*retBox.d<sDistPointToLight) {
				ret.dif *= (1.0 - cub[i].mat.lap)*cub[i].mat.dif;
				ret.spe *= (1.0 - cub[i].mat.lap)*cub[i].mat.spe;
			}
		}
	return ret;
}

material getMaterialFromObjectId(int idObj) {
	if(idObj == floor.id)
		return floor.mat;
	if(idObj == sph.id)
		return sph.mat;
	if(idObj == cub[0].id)
		return cub[0].mat;
	return cub[1].mat;
}

// Implements the Phong model in the given point
vec3 phong(int idObj, vec3 pointObj, vec3 normObj) {
	material mat = getMaterialFromObjectId(idObj);
	// Ambient component
	vec3 color = ambientLight*mat.amb;
	for(int i=0;i<2;++i) {
		vec3 ligDir = normalize(lig[i].pos - pointObj);
		lightComponents ligCom = lightReceived(lig[i], ligDir, idObj, pointObj);
		vec3 color1 = vec3(0.0);
		// Diffuse component
		if(ligCom.dif[0]!=0.0 || ligCom.dif[1]!=0.0 || ligCom.dif[2]!=0.0)
			color1 += mat.dif*max(dot(ligDir, normObj), 0.0)*ligCom.dif;
		// Specular component
		if(ligCom.spe[0]!=0.0 || ligCom.spe[1]!=0.0 || ligCom.spe[2]!=0.0) {
			vec3 ref = reflect(-ligDir, normObj);
			vec3 cam = normalize(cam.pos - pointObj);
			color1 += mat.spe*pow(max(dot(ref, cam), 0.0), mat.shi)*ligCom.spe;
		}
		color += mat.lap*color1;
	}
	return color;
}

float getLightAbsorptionPercentageFromObjectId(int idObj) {
	if(idObj == floor.id)
		return floor.mat.lap;
	if(idObj == sph.id)
		return sph.mat.lap;
	if(idObj == cub[0].id)
		return cub[0].mat.lap;
	return cub[1].mat.lap;
}

float getRefractiveIndexFromObjectId(int idObj) {
	if(idObj == floor.id)
		return floor.mat.rfi;
	if(idObj == sph.id)
		return sph.mat.rfi;
	if(idObj == cub[0].id)
		return cub[0].mat.rfi;
	return cub[1].mat.rfi;
}

// Reflectance computation
float reflectance(float n1, float n2, float cosThetaI, float cosThetaT) {
	float num = n1*cosThetaI - n2*cosThetaT;
	float den = n1*cosThetaI + n2*cosThetaT;
	float rPer = (num*num)/(den*den);
	num = n2*cosThetaI - n1*cosThetaT;
	den = n2*cosThetaI + n1*cosThetaT;
	return 0.5*(rPer + (num*num)/(den*den));
}

// Main ray tracer method
vec3 traceFirstRay(ray r) {
	heapNode heap[HEAPSIZE];
	heap[0].exists = true;
	heap[0].idObj = 0;
	heap[0].r = r;
	// Top-down approach that computes all node colors and rays in the heap
	const int halfSize = HEAPSIZE/2;
	for(int i=0;i<HEAPSIZE;++i)
#define lChild 2*i + 1	// Left child index
#define rChild 2*i + 2	// Right child index
		if(heap[i].exists) {
			detectObjectReturn ret = detectObject(heap[i].r, heap[i].idObj);
			if(ret.inter) {
				heap[i].color = phong(ret.idObj, ret.pointObj, ret.normObj);
				if(i < halfSize) {
					float lap = getLightAbsorptionPercentageFromObjectId(ret.idObj);
					if(lap != 1.0) {
						float n1, n2;
						if(heap[i].idObj != ret.idObj) {
							n1 = AIRRFI;
							n2 = getRefractiveIndexFromObjectId(ret.idObj);
						} else {
							n1 = getRefractiveIndexFromObjectId(heap[i].idObj);
							n2 = AIRRFI;
							ret.normObj = -ret.normObj;
						}
						float cosThetaI = -dot(heap[i].r.d, ret.normObj);
						float eta = n1/n2;
						float sqdSinThetaT = eta*eta*(1.0 - cosThetaI*cosThetaI);
						float refl;
						// Total internal reflection
						if(sqdSinThetaT > 1.0) {
							heap[lChild].exists = false;
							refl = 1.0;
						} else {
							float cosThetaT = sqrt(1.0 - sqdSinThetaT);
							refl = reflectance(n1, n2, cosThetaI, cosThetaT);
							// Refraction
							heap[lChild].exists = true;
							heap[lChild].idObj = ret.idObj;
							heap[lChild].r.o = ret.pointObj;
							heap[lChild].r.d = eta*heap[i].r.d +
								(eta*cosThetaI - cosThetaT)*ret.normObj;
							heap[lChild].mult = (1.0 - lap)*(1.0 - refl);
						}
						// Reflection
						heap[rChild].exists = true;
						heap[rChild].idObj = ret.idObj;
						heap[rChild].r.o = ret.pointObj;
						heap[rChild].r.d = heap[i].r.d + (2.0*cosThetaI)*ret.normObj;
						heap[rChild].mult = (1.0 - lap)*refl;
					} else {
						heap[lChild].exists = false;
						heap[rChild].exists = false;
					}
				}
			} else {
				if(i == 0)
					return BACKGROUND;
				heap[i].color = BACKGROUND;
				if(i < halfSize) {
					heap[lChild].exists = false;
					heap[rChild].exists = false;
				}
			}
		} else
			if(i < halfSize) {
				heap[lChild].exists = false;
				heap[rChild].exists = false;
			}
	// Bottom-up approach that combines the node colors, from the deepests to the first
	for(int i=halfSize;i>=0;--i) {
		if(heap[lChild].exists)
			heap[i].color += heap[lChild].mult*heap[lChild].color;
		if(heap[rChild].exists)
			heap[i].color += heap[rChild].mult*heap[rChild].color;
	}
	return heap[0].color;
}

void main() {
	init();
	vec3 point = cam.pos + cam.lav +
		(gl_FragCoord.x/u_resolution.x - 0.5)*cam.siv +
		(gl_FragCoord.y/u_resolution.y - 0.5)*cam.upv;
	gl_FragColor = vec4(traceFirstRay(ray(cam.pos, normalize(point - cam.pos))), 1.0);
}