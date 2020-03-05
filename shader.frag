// Author:	Santiago Viertel
// Title:	Ray Tracing Shader

precision mediump float;

uniform vec2 u_resolution;

const float EPSILON 	= 0.000001;		// Small constant for avoid overflow
const float AIRRFI		= 1.00027316;	// Refractive index of the air
const vec3 BACKGROUND	= vec3(1.0);	// Background color

struct camera {
	vec3 pos;	// Position
	vec3 upv; 	// Up vector
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

struct sphere {
	int id;			// Object id number
	float rad;		// Radius
	vec3 pos;		// Position
	material mat;	// Material
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

struct ray {
	vec3 o;	// Ray origin
	vec3 d;	// Ray direction
};

struct boxDistanceReturn {
	float d;	// Distance from the origin of the ray
	vec3 n;		// Normal vector of the surface
};

camera cam;
float ambientLight;
pointLight lig[2];
rhombusBox floor;
sphere sph;
rhombusBox cub[2];

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
		1,								// Object id number
		4.0,							// Half of the side of the rhombus
		referenceSystem(
			vec3(-1.8, 1.75751,-1.5),	// Origin of the reference system
			mat3(
				1.0, 0.0, 0.0,			// X vector
				0.0, 1.0, 0.0,			// Y vector
				0.0, 0.0, 1.0			// Z vector (normal)
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
}

boxDistanceReturn boxDistance(ray r, rhombusBox bo) {
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
				if(i==0 || t1>tMin) {
					tMin = t1;
					nMin = bo.rs.v[i];
				}
				if(i==0 || t2<tMax) {
					tMax = t2;
					nMax = -bo.rs.v[i];
				}
			} else {
				if(i==0 || t2>tMin) {
					tMin = t2;
					nMin = -bo.rs.v[i];
				}
				if(i==0 || t1<tMax) {
					tMax = t1;
					nMax = bo.rs.v[i];
				}
			}
			if(tMin>tMax || tMax<=0.0)
				return boxDistanceReturn(-1.0, vec3(0.0));
		} else if(-e-bo.halfSide>0.0 || bo.halfSide-e<0.0)
			return boxDistanceReturn(-1.0, vec3(0.0));
	}
	if(tMin > 0.0)
		return boxDistanceReturn(tMin, nMin);
	return boxDistanceReturn(tMax, nMax);
}

float sphereDistance(ray r, sphere sp) {
	vec3 l = sp.pos - r.o;
	float s1 = dot(l, r.d);
	float ls = dot(l, l);
	float rs = sp.rad*sp.rad;
	float ms = ls - s1*s1;
	if((ls>=rs && s1<=0.0) || ms>=rs)
		//Origin of the ray after the sphere center AND outside the sphere OR
		//Ray above or below the sphere
		return -1.0;
	float q = sqrt(rs - ms);
	if(ls < rs)
		// Origin of the ray inside the sphere
		return s1 + q;
	return s1 - q;
}

// http://graphics.cs.kuleuven.be/publications/LD04ERQIT/LD04ERQIT_paper.pdf
float rhombusDistance(ray r, rhombusBox ro) {
	vec3 v00 = ro.rs.o + ro.halfSide*(ro.rs.v[0] + ro.rs.v[1]);
	vec3 v10 = ro.rs.o + ro.halfSide*(ro.rs.v[0] - ro.rs.v[1]);
	vec3 v01 = ro.rs.o + ro.halfSide*(ro.rs.v[1] - ro.rs.v[0]);
	vec3 e01 = v10 - v00;
	vec3 e03 = v01 - v00;
	vec3 vect1 = cross(r.d, e03);
	float det = dot(e01, vect1);
	if(abs(det) < EPSILON)
		return -1.0;
	vec3 vect2 = r.o - v00;
	float alpha = dot(vect2, vect1)/det;
	if(alpha<0.0 || alpha>1.0)
		return -1.0;
	vec3 vect3 = cross(vect2, e01);
	float beta = dot(r.d, vect3)/det;
	if(beta<0.0 || beta>1.0)
		return -1.0;
	float t = dot(e03, vect3)/det;
	if(t <= 0.0)
		return -1.0;
	return t;
}

lightComponents lightReceived(pointLight lig, vec3 ligDir, int objId, vec3 objPoint) {
	ray r = ray(objPoint, ligDir);
	lightComponents ret = lig.com;
	vec3 pointToLight = lig.pos - objPoint;
	float sDistFromPointToLight = dot(pointToLight, pointToLight);
	float d;
	if(objId != floor.id) {
		d = rhombusDistance(r, floor);
		if(d!=-1.0 && d*d<sDistFromPointToLight) {
			ret.dif *= (1.0 - floor.mat.lap)*floor.mat.dif;
			ret.spe *= (1.0 - floor.mat.lap)*floor.mat.spe;
		}
	}
	if(objId != sph.id) {
		d = sphereDistance(r, sph);
		if(d!=-1.0 && d*d<sDistFromPointToLight) {
			ret.dif *= (1.0 - sph.mat.lap)*sph.mat.dif;
			ret.spe *= (1.0 - sph.mat.lap)*sph.mat.spe;
		}
	}
	for(int i=0;i<2;++i)
		if(objId != cub[i].id) {
			boxDistanceReturn retCub = boxDistance(r, cub[i]);
			if(retCub.d!=-1.0 && retCub.d*retCub.d<sDistFromPointToLight) {
				ret.dif *= (1.0 - cub[i].mat.lap)*cub[i].mat.dif;
				ret.spe *= (1.0 - cub[i].mat.lap)*cub[i].mat.spe;
			}
		}
	return ret;
}

material getMaterialFromObjectId(int objId) {
	if(objId == floor.id)
		return floor.mat;
	if(objId == sph.id)
		return sph.mat;
	if(objId == cub[0].id)
		return cub[0].mat;
	return cub[1].mat;
}

vec3 shade(vec3 camDir, int objId, vec3 objPoint, vec3 objNorm) {
	material mat = getMaterialFromObjectId(objId);
	// Ambient component
	vec3 color = ambientLight*mat.amb;
	for(int i=0;i<2;++i) {
		vec3 ligDir = normalize(lig[i].pos - objPoint);
		lightComponents ligCom = lightReceived(lig[i], ligDir, objId, objPoint);
		vec3 color1 = vec3(0.0);
		// Diffuse component
		if(ligCom.dif[0]!=0.0 || ligCom.dif[1]!=0.0 || ligCom.dif[2]!=0.0)
			color1 += mat.dif*max(dot(ligDir, objNorm), 0.0)*ligCom.dif;
		// Specular component
		if(ligCom.spe[0]!=0.0 || ligCom.spe[1]!=0.0 || ligCom.spe[2]!=0.0) {
			vec3 ref = reflect(-ligDir, objNorm);
			color1 += mat.spe*pow(max(dot(ref, camDir), 0.0), mat.shi)*ligCom.spe;
		}
		color += mat.lap*color1;
	}
	return color;
}

vec3 traceRay(ray r) {
	int objId = 0;
	// Floor
	float objDist = rhombusDistance(r, floor);
	if(objDist != -1.0)
		objId = floor.id;
	// Sphere
	float d = sphereDistance(r, sph);
	if(d!=-1.0 && (objDist==-1.0 || d<objDist)) {
		objId = sph.id;
		objDist = d;
	}
	vec3 objNorm;
	// Boxes
	for(int i=0;i<2;++i) {
		boxDistanceReturn ret = boxDistance(r, cub[i]);
		if(ret.d!=-1.0 && (objDist==-1.0 || ret.d<objDist)) {
			objId = cub[i].id;
			objDist = ret.d;
			objNorm = ret.n;
		}
	}
	if(objId == 0)
		return BACKGROUND;
	vec3 objPoint = r.o + objDist*r.d;
	if(objId == floor.id)
		objNorm = floor.rs.v[2];
	else if(objId == sph.id)
		objNorm = normalize(objPoint - sph.pos);
	return shade(-r.d, objId, objPoint, objNorm);
}

void main() {
	init();
	vec3 ponto = cam.pos + cam.lav + (gl_FragCoord.x/u_resolution.x - 0.5)*cam.siv + (gl_FragCoord.y/u_resolution.y - 0.5)*cam.upv;
	gl_FragColor = vec4(traceRay(ray(cam.pos, normalize(ponto - cam.pos))), 1.0);
}