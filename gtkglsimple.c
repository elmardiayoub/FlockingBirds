#include <stdio.h>
#include <stdlib.h>
#include <gtk/gtk.h>
#include <gtkgl/gtkglarea.h>
#include <GL/gl.h>
#include <GL/glu.h>

#define DISTANCE_INIT 6.0
#define DISTANCE_MIN 4.5
#define DISTANCE_MAX 25.0

////////////////////////////////////////////////////////////////////////////////////////////////////////



#include <stdbool.h>


#include <sys/time.h>
#include <unistd.h>
#include<math.h>




#define PI 3.14159265



///////////////////////////////////////////////////////////////////////////////////////////

#define NBbird 2000
#define NBbirds 2000
#define NBpoint 2000
#define NBvector 2000
#define NBtile 2000

#define TSize 1.0
#define GridSize 40

typedef struct liste
{
  double x;
  double y;
  struct liste *svt;
}lis;

typedef struct ps
{
  int ind;
  struct ps *svt;
}speedFile;


typedef struct 
{
  lis* tete;
  lis* queue;

}file;

typedef struct 
{
  speedFile* tete;
  speedFile* queue;

}file1;

int tmp=0,tmp2=0;
int id=0;
file* mafile=NULL;
file1 *File2 =NULL;
int play=0;
int start=0;
int nbspeed=0;

enum direction { TOP , TOP_RIGHT, RIGHT, BOTTOM_RIGHT, BOTTOM, BOTTOM_LEFT, LEFT, TOP_LEFT,CURRENT};
/// MAP :
typedef struct
{
    int tab[NBbird];
    int IDE;

}Mybirdtab;

typedef struct tile
{
    int tileIndex;
    bool ocupated;
    Mybirdtab* localBirds;
    int nearObstacle;
}MyTile;

typedef struct
{
    MyTile* tab[NBtile];
    int IDE;

}Mytiletab;


typedef struct maps
{
    float wsPosX;	///-20.0
    float wsPosZ;   ///-20.0
    int xTiles;
    int yTiles;
    float tileSize;

    Mytiletab* tiles;

}MyMap;


MyTile* init_tile();
void inserer_tab(Mytiletab* tableau, MyTile* tile);
MyMap* init_map(int xTiles,int yTiles,float tileSize);

MyTile* getTileByGridCoords(MyMap* mapp,double x, double y);
void draw_map(MyMap* mapp);


///Tpoint :

typedef struct
{
    double x,y,z;
}TPoint;

    TPoint* init_point(double x, double y, double z);
    TPoint* add(TPoint* p1,TPoint* p2);
    TPoint* inverser(TPoint* p);
    TPoint* lerp(TPoint* p,TPoint* p2,float timeStep);


///Tvector  :
typedef struct
{
    double x,y,z;
    TPoint* origin;

}Tvector;

    Tvector* init_vector(double x, double y, double z);
    Tvector* init_vector_with_two_points(TPoint* p1, TPoint* p2);
    Tvector* init_vector_with_one_points(TPoint* p1);
    double magnitude(Tvector* v);
    Tvector* unit(Tvector* v);
    Tvector* add_vector(Tvector* v1, Tvector* v2);
    Tvector* subtract(Tvector* v1, Tvector* v2);
    double dotProduct(Tvector* v1, Tvector* v2);
    double dotProduct_angle(Tvector* v1, Tvector* v2,double angle);
    double getAngle(Tvector* v1,Tvector* v2);
    Tvector* crossProductP(Tvector* v1, Tvector* v2);
    Tvector* crossProduct(Tvector* v1, Tvector* v2);
    Tvector* multiply(Tvector* v1, float s);
    double  scalarTripleProduct(Tvector* v1, Tvector* v2, Tvector* v3);
    Tvector* normalCalc(Tvector* v1,Tvector* v2,Tvector* v3);
    void draw_vector(Tvector* v);




///Bird :

typedef struct
{
    TPoint* tab[NBpoint];
    int IDE;

}Tpointab;

typedef struct
{
    Tvector tab[NBvector];
    int IDE;

}Tvectortab;



typedef struct
{

  TPoint* pos;
	Tvector* speed;
	int tileIndex;
	float angle;
	bool inColision;
	double age;
  Tpointab* coh;
  Tvectortab* algn;
  Tvectortab* sep;
  Tpointab* nextPos;
  int ColorIset;
  double red;
  double green;
  double blue;

}Bird;



typedef struct
{
    Bird* tab[NBbirds];
    int IDE;

}birdtab;


void draw_bird(Bird* bird);
void inserer_tab_bird(birdtab* tableau, Bird* bird);
Tvector* steer;
Tvector* vect1;
Tvector* vect2;
Tvector* vect3;
Tvector* vect4;
Tvector* vect5;
Tvector* vect6;
Tvector* vect7;
Tvector* vect8;
Tvector* vect9;
Tvector* vect10;
Tvector* vect11;
Tvector* vect12;
Tvector* sepTmp;
TPoint* cohTmp;
Tvector* algnTmp;
TPoint *Point;


/////////////////////////////////////////////////////////////////////////////////////////////////////


float power=0.02;
float force=0.24;
MyTile* init_tile()
{
    MyTile* tile;

    tile=(MyTile*)malloc(sizeof(MyTile));
    if(!tile)
    {
        printf("Erreur 1\n");
        exit(-1);
    }

    tile->tileIndex=0;
    tile->ocupated=FALSE;
    tile->nearObstacle=-1;

    tile->localBirds=(Mybirdtab*)malloc(sizeof(Mybirdtab));
    if(!tile->localBirds)
    {
        printf("Erreur 2\n");
        exit(-1);
    }
    tile->localBirds->IDE=-1;
    return((MyTile*)tile);
}

void inserer_tab(Mytiletab* tableau, MyTile* tile)
{


    tableau->tab[++tableau->IDE]=tile;

}

void occupy(MyMap* map,int index)
{
	map->tiles->tab[index]->ocupated=TRUE;
}

void inserer_localbird(Mybirdtab* tableau, int birdID)
{

    tableau->tab[++tableau->IDE]=birdID;

}

void addBoid(MyMap* map,int index,int birdId)
{


    inserer_localbird(map->tiles->tab[index]->localBirds,birdId);



}

int getXGridCoord(MyMap* map,int index)
{
	
	int xCoord;
	xCoord = index%map->xTiles;
	return xCoord;
}

int getYGridCoord(MyMap* map,int index)
{
	int yCoord;
	yCoord= index/map->xTiles;

	return yCoord;
}
MyTile* getTileByIndex(MyMap* map,int i)
{
	return map->tiles->tab[i];
}

MyTile* getTileByGridCoords(MyMap* map,double x, double y)
{
	int index = (y * map->xTiles) + x;
	return map->tiles->tab[index];
}

MyTile *getNeighbouringTile(MyMap* map,int index, int direction)
{
	int xPos = getXGridCoord(map,index);
	int yPos = getYGridCoord(map,index);
	MyTile* tile = getTileByIndex(map,index);

	switch (direction)
	{
		case TOP:
		
			if (yPos == 0 ) return NULL;	
			return getTileByGridCoords(map,xPos, yPos - 1);
			break;

		case TOP_RIGHT:
			if (yPos ==0 || xPos == map->xTiles-1 ) return NULL;
			return getTileByGridCoords(map,xPos+1, yPos - 1);
			break;

		case RIGHT:
			if (xPos == map->xTiles-1 ) return NULL;
			return getTileByGridCoords(map,xPos+1, yPos );
			break;

		case BOTTOM_RIGHT:
			if (yPos == map->yTiles-1 || xPos == map->xTiles-1 ) return NULL;
			return getTileByGridCoords(map,xPos+1, yPos + 1);
			break;

		case BOTTOM:
			if (yPos == map->yTiles-1 ) return NULL;
			return getTileByGridCoords(map,xPos, yPos + 1);
			break;

		case BOTTOM_LEFT:
			if (yPos == map->yTiles-1 || xPos == 0 ) return NULL;
			return getTileByGridCoords(map,xPos-1, yPos + 1);
			break;

		case LEFT:
			if (xPos == 0 ) return NULL;
			return getTileByGridCoords(map,xPos-1, yPos);
			break;

		case TOP_LEFT:
			if (yPos == 0 || xPos == 0 ) return NULL;
			return getTileByGridCoords(map,xPos-1, yPos - 1);
			break;
		case CURRENT:
			return getTileByIndex(map,index);
			break;
	}
	return NULL;
}


MyMap* init_map(int xTiles,int yTiles,float tileSize)
{
    MyMap* map;
    MyTile* t = NULL;

    int x,y;


    map=(MyMap*)malloc(sizeof(MyMap));
    if(!map)
    {
        printf("Erreur 3\n");
        exit(-1);
    }
    map->tiles = (Mytiletab*)malloc(sizeof(Mytiletab));
    if(!map->tiles)
    {
        printf("Erreur 4\n");
        exit(-1);
    }
    map->tiles->IDE=-1;
    map->xTiles = xTiles;
    map->yTiles = yTiles;
    map->tileSize = tileSize;


    for(y=0;y<yTiles;y++)
    {
        for(x=0;x<xTiles;x++)
        {
            t=init_tile();
            t->tileIndex = x + y * xTiles;
            inserer_tab(map->tiles,t);
        }
    }

    map->wsPosX=-20.0;
    map->wsPosZ=-20.0;

    return((MyMap*)map);

}


void draw_map(MyMap* mapp)
{
   	MyTile *t = NULL;


	glTranslatef(mapp->wsPosX, 0.0, mapp->wsPosZ);
	for (int z = 0; z < mapp->yTiles; z++)
	{
		for (int x = 0; x < mapp->xTiles; x++)
		{
			t = getTileByGridCoords(mapp,x,z);

            glColor3f(0.573, 0.573, 0.573);
            glBegin(GL_LINE_STRIP);
                glVertex3f((float) x * mapp->tileSize, 0.0, (float) z * mapp->tileSize);
                glVertex3f((float) x * mapp->tileSize, 0.0, (float)(z * mapp->tileSize) + mapp->tileSize);
                glVertex3f((float)(x * mapp->tileSize) + mapp->tileSize, 0.0, (float)(z * mapp->tileSize) + mapp->tileSize);
                glVertex3f((float)(x * mapp->tileSize) + mapp->tileSize, 0.0, (float) z * mapp->tileSize);
                glVertex3f((float) x * mapp->tileSize, 0.0, (float) z * mapp->tileSize);
            glEnd();
		}
	}
}

///****************** - TPoint  - *********************///:


TPoint* init_point(double x, double y, double z)
{
    TPoint* p;

    p=(TPoint*)malloc(sizeof(TPoint));
    if(!p)
    {
        printf("Erreur 5\n");
        exit(-1);
    }
    p->x=x;
    p->y=y;
    p->z=z;

    return((TPoint*)p);
}
TPoint* init_point_2(TPoint* p,double x, double y, double z)
{

    p->x=x;
    p->y=y;
    p->z=z;

    return((TPoint*)p);
}

void inserer_tab_point(Tpointab* tableau, TPoint* point)
{

    tableau->tab[++tableau->IDE]=point;

}


TPoint* add(TPoint* p1,TPoint* p2)
{

    return((TPoint*)init_point((p1->x)+(p2->x),(p1->y)+(p2->y),(p1->z)+(p2->z)));
}
TPoint* add_2(TPoint* p1,TPoint* p2)
{
    p1->x=(p1->x)+(p2->x);
    p1->y=(p1->y)+(p2->y);
    p1->z=(p1->z)+(p2->z);
    return((TPoint*)p1);
}
TPoint* add_point_1p(TPoint* p,TPoint* p1,TPoint* p2)
{
    p->x=(p1->x)+(p2->x);
    p->y=(p1->y)+(p2->y);
    p->z=(p1->z)+(p2->z);
    return((TPoint*)p);
}
TPoint* add_pos(TPoint* p1,TPoint* p2)
{
    Point->x=(p1->x)+(p2->x);
    Point->y=(p1->y)+(p2->y);
    Point->z= (p1->z)+(p2->z);
    return((TPoint*)Point);
}
TPoint* inverser_pos(TPoint* p)
{
    Point->x = -p->x;
    Point->y = -p->y;
    Point->z = -p->z;
    return((TPoint*)Point);
}
TPoint* inverser(TPoint* p)
{
    return((TPoint*)init_point(-p->x,-p->y,-p->z));
}

TPoint* lerp(TPoint* p,TPoint* p2,float timeStep)
{
	TPoint* p1=init_point(p->x,p->y,p->z);
	TPoint* delta=add(p2,inverser(p1));
	delta->x=delta->x*timeStep;
	delta->y=delta->y*timeStep;
	delta->z=delta->z*timeStep;

	return add(p1,delta);
}


///****************** - TVector  - *********************///:

Tvector* init_vector(double x, double y, double z)
{
    Tvector* v;

    v=(Tvector*)malloc(sizeof(Tvector));
    if(!v)
    {
        printf("Erreur 6\n");
        exit(-1);
    }
    v->x=x;
    v->y=y;
    v->z=z;

    v->origin=init_point(0,0,0);

    return((Tvector*)v);
}

void inserer_tab_vector(Tvectortab* tableau, Tvector* vec)
{

    tableau->tab[++tableau->IDE]=*vec;

}

Tvector* init_vector_with_two_points(TPoint* p1, TPoint* p2)
{
    Tvector* v;

    v=(Tvector*)malloc(sizeof(Tvector));
    if(!v)
    {
        printf("Erreur 7\n");
        exit(-1);
    }
    v->x=p1->x;
    v->y=p1->y;
    v->z=p1->z;

    v->origin=p2;

    return((Tvector*)v);
}
Tvector* init_vector_with_one_points(TPoint* p1)
{
    Tvector *v;

    v=(Tvector*)malloc(sizeof(Tvector));
    if(!v)
    {
        printf("Erreur 8\n");
        exit(-1);
    }
    v->x=p1->x;
    v->y=p1->y;
    v->z=p1->z;

    v->origin=init_point(0,0,0);

    return((Tvector*)v);
}
Tvector* init_vector_1p(Tvector* v1,TPoint* p1)
{

    v1->x=p1->x;
    v1->y=p1->y;
    v1->z=p1->z;

    return((Tvector*)v1);
}
Tvector* init_vector_steer(TPoint* p1)
{

    steer->x=p1->x;
    steer->y=p1->y;
    steer->z=p1->z;

    return((Tvector*)steer);
}

double magnitude(Tvector* v)
{
	return (double)(sqrt(pow(v->x,2)+pow(v->y,2)+pow(v->z,2))); //directly return the calculated magnitude without storing it
}


Tvector* unit(Tvector* v)
{

	double x,y,z;
	double mag;


	mag=(double)(sqrt(pow(v->x,2)+pow(v->y,2)+pow(v->z,2)));

	x=(double)v->x/mag; 
	y=(double)v->y/mag;
	z=(double)v->z/mag;
	return ((Tvector*)init_vector(x,y,z)); 
}
Tvector* unit_2(Tvector* v)
{

	double mag;


	mag=(double)(sqrt(pow(v->x,2)+pow(v->y,2)+pow(v->z,2)));

   if(mag)
   {
      v->x=(double)v->x/mag; 
      v->y=(double)v->y/mag;
      v->z=(double)v->z/mag;
   }

	return ((Tvector*)v); //return a TVector type item
}
Tvector* unit_vect(Tvector* v1,Tvector* v2)
{

	double mag;


	mag=(double)(sqrt(pow(v2->x,2)+pow(v2->y,2)+pow(v2->z,2)));

	 // get the vector magnitude
	v1->x=(double)v2->x/mag; // normalize each component
	v1->y=(double)v2->y/mag;
	v1->z=(double)v2->z/mag;
	return ((Tvector*)v1); //return a TVector type item
}


Tvector* add_vector(Tvector* v1, Tvector* v2)
{

	//add each component from the 2 vectors
	v1->x=v1->x+v2->x;
	v1->y=v1->y+v2->y;
	v1->z=v1->z+v2->z;
	return ((Tvector*)v1); //return TVector type item
}


Tvector* subtract(Tvector* v1, Tvector* v2)
{
	double x,y,z;
	//add each component from the 2 vectors
	x=v1->x-v2->x;
	y=v1->y-v2->y;
	z=v1->z-v2->z;
	return ((Tvector*)init_vector(x,y,z)); //return Tvector type item
}
Tvector* subtract_2(Tvector* v1, Tvector* v2)
{
//	double x,y,z;
//	//add each component from the 2 vectors
	v1->x=v1->x-v2->x;
	v1->y=v1->y-v2->y;
	v1->z=v1->z-v2->z;
  
	return ((Tvector*)v1); //return Tvector type item
}

double dotProduct(Tvector* v1, Tvector* v2)
{
	double result;
	result=(v1->x)*(v2->x)+(v1->y)*(v2->y)+(v1->z)*(v2->z); //calculate result and the return
	return result;
}
double dotProduct_angle(Tvector* v1, Tvector* v2,double angle)
{
	angle=angle*PI/180; // Transform the angle in radians
	return (magnitude(v1)*magnitude(v2)*cos(angle));
}

double getAngle(Tvector* v1,Tvector* v2)
{

	double angle;
	angle=acos(dotProduct(unit_2(v1),unit_2(v2)));
	    //calculate the angle in radians
	angle=(double)angle*180/PI;

    // transform the angle in degrees
	return (double)angle;
}

Tvector* crossProductP(Tvector* v1, Tvector* v2)
{
	double x,y,z;
	//calculate each component
	x=abs(v1->y)*abs(v2->z)+abs(v1->z)*abs(v2->y);
	y=abs(v1->z)*abs(v2->x)+abs(v1->x)*abs(v2->z);
	z=abs(v1->x)*abs(v2->y)+abs(v1->y)*abs(v2->x);

	return ((Tvector*)init_vector(x,y,z));
}


Tvector* crossProduct(Tvector* v1, Tvector* v2)
{
//	double x,y,z;
	//calculate each component
	v1->x=v1->y*v2->z-v1->z*v2->y;
	v1->y=v1->z*v2->x-v1->x*v2->z;
	v1->z=v1->x*v2->y-v1->y*v2->x;

	return ((Tvector*)v1);
}

Tvector* multiply(Tvector* v1, float s)
{
    v1->x=(double)(v1->x)*s;
    v1->y=(double)(v1->y)*s;
    v1->z=(double)(v1->z)*s;

	return ((Tvector*)v1);
}

double  scalarTripleProduct(Tvector* v1, Tvector* v2, Tvector* v3)
{
	double stp;
	//calculate the dot product between a vector and the cross product of the other 2 vectors
	stp=dotProduct(v1,crossProduct(v2,v3));
	return stp;
}
Tvector* normalCalc(Tvector* v1,Tvector* v2,Tvector* v3)
{

    return ((Tvector*)unit(crossProduct(subtract(v1,v3),subtract(v2,v3))));
}
void draw_vector(Tvector* v)
{
	glBegin(GL_LINES);
	glVertex3f( v->origin->x, v->origin->y, v->origin->z);
    glVertex3f(v->x, v->y, v->z);
    glEnd();
}

///****************** - Bird  - *********************///:




void inserer_tab_bird(birdtab* tableau, Bird* bird)
{

    tableau->tab[++tableau->IDE]=bird;

}

Bird* init_bird(TPoint* pos, Tvector* speed,double age)
{
    Bird* b;
    b=(Bird*)malloc(sizeof(Bird));
    if(!b)
    {
        printf("Erreur 9\n");
        exit(-1);
    }
    b->algn=(Tvectortab*)malloc(sizeof(Tvectortab));
    if(!b->algn)
    {
        printf("Erreur 10\n");
        exit(-1);
    }
    b->coh=(Tpointab*)malloc(sizeof(Tpointab));
    if(!b->coh)
    {
        printf("Erreur 11\n");
        exit(-1);
    }
    b->sep=(Tvectortab*)malloc(sizeof(Tvectortab));
    if(!b->sep)
    {
        printf("Erreur 12\n");
        exit(-1);
    }
    b->nextPos=(Tpointab*)malloc(sizeof(Tpointab));
    if(!b->nextPos)
    {
        printf("Erreur 13\n");
        exit(-1);
    }
    b->ColorIset=0;
    b->pos=pos;

	b->speed=unit(speed);
//    printf("[%lf]\n",b->speed);


	//this->speed.unit();
	int x,z;
	x=(int)(pos->x/2);
	z=(int)(pos->z/2);
	b->tileIndex=z*20+x;
	b->angle=getAngle(init_vector(0,0,1),speed);
	if(crossProduct(init_vector(0,0,1),speed)->y<0) b->angle=360-b->angle;
	b->inColision=FALSE;
	b->age=age;

    b->sep->IDE=-1;
    b->algn->IDE=-1;
    b->coh->IDE=-1;
    b->nextPos->IDE=-1;

	return((Bird*)b);
}

void clear_bird(Bird* b)
{
    b->sep->IDE=-1;
    b->algn->IDE=-1;
    b->coh->IDE=-1;
    b->nextPos->IDE=-1;
}



void draw_bird(Bird* bird)
{

  if(bird->ColorIset)
    glColor3f(bird->red,bird->green,bird->blue);
  
	glPushMatrix();
	//glColor3f(0.0,1.0,0.0);
	//if(this->inColision==true) glColor3f(1.0,0,0);
  // printf("%lf   %lf    %lf  \n",bird->pos->x,bird->pos->y,bird->pos->z);
  	glTranslatef(-20.0, 0.0, -20.0);
	glTranslatef(bird->pos->x,bird->pos->y,bird->pos->z);
	glRotatef(bird->angle,0.0,1,0);
	glBegin(GL_TRIANGLES);		// Drawing Using Triangles
	glVertex3f( 0.0, 0, 0.5*bird->age);		// Top
	glVertex3f(-0.2*bird->age, 0, -0.3*bird->age);		// Bottom Left
	glVertex3f( 0.2*bird->age, 0, -0.3*bird->age);		// Bottom Right
    glEnd();
	glPopMatrix();
	//TVector(pos,TPoint(pos._x+speed._x,pos._y+speed._y,pos._z+speed._z)).draw();

}

int i=0;
///
void addPosition(Bird* b,TPoint* p)
{
    inserer_tab_point(b->nextPos, p);
}
Tvector* init(Tvector* v)
{
    v->x=0;
    v->y=0;
    v->z=1;
    return((Tvector*)v);
}
Tvector* init_2(Tvector* v)
{
    v->x=0;
    v->y=0;
    v->z=0;
    return((Tvector*)v);
}

TPoint* init_2p(TPoint* p)
{
    p->x=0;
    p->y=0;
    p->z=0;
    return((TPoint*)p);
}

void updateSpeed(Bird* b,Tvector* avoidance)
{
//    printf("Avant : %lf",b->speed->x);
//    i++;
//    printf("%f\n",force);
//       printf(" %lf\n",avoidance->z);
//    printf("%lf ---- \n",multiply(avoidance,(float)power)->z);
//        printf("%f\n",power);
	b->speed=add_vector(b->speed,multiply(avoidance,(float)power));
//	printf("%lf ---- \n",multiply(avoidance,(float)power)->x);

    double mag;

//	printf("%lf\n",v->x);
	mag=(double)(sqrt(pow(b->speed->x,2)+pow(b->speed->y,2)+pow(b->speed->z,2)));
//	printf("%lf\n",mag);
	 // get the vector magnitude

//	  printf("%lf \n",pow(b->speed->z,2));

	b->speed->x=(double)(b->speed->x)/(double)mag; // normalize each component
	b->speed->y=(double)(b->speed->y)/(double)mag;
	b->speed->z=(double)(b->speed->z)/(double)mag;

//	return ((Tvector*)init_vector(x,y,z));

//    printf("%lf \n ",b->speed->x);
//   printf(" %lf \n ",b->speed->x);
//	b->speed=unit(b->speed);
//	printf(" ------ %lf \n",b->speed->x);

	b->angle=getAngle(init(vect9),b->speed);
//	printf("%lf ---\n ",vect3->x);
//	 printf("%lf ---\n ",b->angle);
	if((crossProduct(init(vect10),b->speed)->y)<0) b->angle=360-(b->angle);
}
double distance(TPoint* p1,TPoint* p2)
{
    double result;
    result=sqrt(pow(p2->x-p1->x,2)+pow(p2->y-p1->y,2)+pow(p2->y-p1->y,2));
//    printf("%lf\n",result);
    return result;
}

void setTileIndex(Bird* b)
{
	int x,z;
	x=b->pos->x/TSize;
	z=b->pos->z/TSize;
  ///une formule pour avoir les indices des carreau
  ///successives
	b->tileIndex=z*GridSize+x;
}
void setNewPos(Bird* b)
{


  //Factor c’est une variable pour diminuer la vitesse
	float factor=0.08;
	float maxLimit=39.5;
	float minLimit=0.5;


	b->pos->x=b->pos->x+(b->speed->x)*((float)factor);
	b->pos->z=b->pos->z+(b->speed->z)*((float)factor);

  ///S’il a dépassé la taille de la fenêtre
	if((b->pos->x)>(maxLimit)) b->pos->x=minLimit;
	else if((b->pos->x)<(minLimit)) b->pos->x=maxLimit;

	if((b->pos->z)>(maxLimit)) b->pos->z=minLimit;
	else if((b->pos->z)<(minLimit)) b->pos->z=maxLimit;

  ///Après le changement de position on change l’indice du tile
  ///ou l’oiseau existe
	setTileIndex(b);

}

void separation(Bird *b,Bird *b1)
{
//    printf("%lf\n",b1->sep);
//    for(int l=0;l<b->sep->IDE+1;l++) printf("%d %p",l,b->sep->tab[i]); printf("\n");
    inserer_tab_vector(b->sep,unit_2(init_vector_1p(vect6,add_point_1p(Point,b->pos,inverser(b1->pos)))));
//

}



void alignment(Bird *b,Bird *b1)
{
	inserer_tab_vector(b->algn,unit_vect(vect5,b1->speed));
}

void cohesion(Bird *b,Bird *b1)
{

	 inserer_tab_point(b->coh,b1->pos);
//	  printf("IDE :  %d\n",b->coh->IDE);
}

void steerTowards(Bird* b,TPoint* futureLoc)
{

//          printf("%lf --  ",steer->x);
//        printf("%lf ---- ",steer->z);
		steer=init_vector_steer(add_pos(inverser_pos(b->pos),futureLoc));
		// printf("%lf \n",steer->x);

//		 printf("%lf \n ",steer->x);
//		printf("%p\n",steer);
//


//
  //  printf("%lf \n",steer->x); 
		steer=subtract_2(unit_2(steer),unit_2(b->speed));
     
//    steer->x=unit(steer)->x-unit(b->speed)->x;
//	steer->y=unit(steer)->y-unit(b->speed)->y;
//
//	steer->z=unit(steer)->z-unit(b->speed)->z;



//    printf("%f\n",power);
        power=(float)0.02;
		updateSpeed(b,steer);

		// printf(" %lf \n",steer->x);


}

void newSpeed(Bird* b,bool useIt)
{

//	Tvector* newSpeed;

    sepTmp=init_2(sepTmp);
    cohTmp=init_2p(cohTmp);
    algnTmp=init_2(algnTmp);

	if(b->coh->IDE+1>0)
	{
		float count=b->sep->IDE+1;
		//nombredes élements de la liste
		for(int i=0;i<b->sep->IDE+1;i++)
		{
//		    printf("%d %lf",i,b->sep->tab[i]->x);
			sepTmp=add_vector(sepTmp,&b->sep->tab[i]);
			cohTmp=add_2(cohTmp,b->coh->tab[i]);
			algnTmp=add_vector(algnTmp,&b->algn->tab[i]);
		}
//		printf("\n");
//        printf("%lf : %lf : %lf \n",sepTmp->x,cohTmp->x,algnTmp->x);
		sepTmp=multiply(sepTmp,(float)1.0/count);

		count=b->algn->IDE+1;

		algnTmp=unit_2(algnTmp);

		Tvector* cohVec=init_vector_1p(vect7,cohTmp);
		cohVec=multiply(cohVec,(float)1.0/count);

		cohVec=subtract_2(cohVec,init_vector_1p(vect11,b->pos));

		if(useIt)
		{
		    power=(float)0.24;
			updateSpeed(b,sepTmp);
			power=(float)0.1;
			updateSpeed(b,algnTmp);
			power=(float)0.165;
			updateSpeed(b,cohVec);
			b->angle=getAngle(init(vect3),b->speed);
			if((crossProduct(init(vect8),b->speed)->y)<0) b->angle=360-b->angle;

    }

	}

	clear_bird(b);
}







//////////////////////////////////////////////////////////////////////////////////////////////////////



birdtab* flock=NULL;
MyMap* map=NULL;

Mybirdtab* usedTiles;
TPoint *averagePos=NULL;

int ind=0;






/*Cette fonction représente le cœur du déplacement des oiseaux car
dans cette fonction on change les positions des oiseaux ainsi que
leur angle, leur séparation, alignement, cohésion et la vitesse.*/
void flockMovement()
{

  usedTiles->IDE=-1;
//  map->tiles->IDE=-1;
 for(int i=0;i<map->tiles->IDE+1;i++)
    map->tiles->tab[i]->localBirds->IDE=-1;

///Pour l'ensemble des oiseaux qu'ils existent
///on change la position et on ajoute au carreau du carte
///l'indice de l’oiseau qu'il l'occupe
  for(int i=0;i<flock->IDE+1;i++)
  {

	  setNewPos(flock->tab[i]);
//    printf("Avt : %d --> %d --> %lf ----\n",ind,i,flock->tab[i]->pos->x);
	  int boidTile=flock->tab[i]->tileIndex;


//     printf("%d %d %d\n",ind,i,usedTiles->tab[i]);
	  addBoid(map,boidTile,i);


  }


///C'est le point du centre de l'ensemble des oiseaux auquels ils se rassemblent
  averagePos=init_point_2(averagePos,0,0,0);
//
  for(int i=0;i<flock->IDE+1;i++)
  {
    ///J'accumule l'ensemble des points des oiseaux
	  averagePos = add_2(averagePos,flock->tab[i]->pos);

///Pour chaque oiseau je parcours les 8 cotés qu'il a (Buttom,top,left,right ...)
	  MyTile *tile;

	  for(int j=0;j<9;j++)
	  {

          //pour chaque coté il me retourne le carreau qu'il correspond
          tile=getNeighbouringTile(map,flock->tab[i]->tileIndex,j);

          ///si c'était pas une cas limite
          ///je pacours l'ensemble des voisins qui existent dans ce carreau
          if (tile!=NULL)
          for(int k=0;k<tile->localBirds->IDE+1;k++)
          {

          /// Je calcule la distance entre le oiseau courant et ces voisins ..
          float result=sqrt(pow(flock->tab[tile->localBirds->tab[k]]->pos->x-flock->tab[i]->pos->x,2));

              //si l'élement courant est différent du voisin c'est à dire lui meme
              //et si angle entre lui et son voisin inférieur à 135
              ///et si la distance entre lui et son voisin inférieur à 1
              ///donc on applique sur ce oiseau les 3 règles
              if(i!=tile->localBirds->tab[k] &&
                getAngle(flock->tab[i]->speed,subtract_2(init_vector_1p(vect1,flock->tab[tile->localBirds->tab[k]]->pos),init_vector_1p(vect2,flock->tab[i]->pos)))<=135 &&
                result<1)
              {
                  //Enregistrer les vecteurs de séparation
                  separation(flock->tab[i],flock->tab[tile->localBirds->tab[k]]);



                  //Enregistrer les vecteurs de l'alignement
                  alignment(flock->tab[i],flock->tab[tile->localBirds->tab[k]]);


                  //Enregistrer les positions des voisins
                  cohesion(flock->tab[i],flock->tab[tile->localBirds->tab[k]]);

              }
		    }
	  }
      ///Changer la vitesse selon les nouvelles changement appliqué sur le tableau de séparation
      ///alignemenet et cohésion
      newSpeed(flock->tab[i],TRUE);

  }

  ///après qu'on _ accumulé tous les points on divise sur le
  //nombre des oiseaux pour avoir
  //une point centre
  averagePos->x=(averagePos->x/(double)(flock->IDE+1));

  averagePos->z=(averagePos->z/(double)(flock->IDE+1));


  ///SteerTowards c'est une fonction qui va permettre
  //aux oiseaux de se diriger au point centre qu'on a
  //calculer
  for(int ind1=0;ind1<flock->IDE+1;ind1++)
       if(flock->tab[ind1]->coh->IDE+1<=1)
            steerTowards(flock->tab[ind1],averagePos);

}



void init_drawing(void)
{

  glShadeModel(GL_FLAT);

  glDisable(GL_LIGHTING);

  glDepthFunc(GL_LEQUAL);
  glEnable(GL_DEPTH_TEST);


  mafile = (file*)malloc(sizeof(file));
  if(!mafile)
  {
    printf("ERREUR 1");
    exit(-1);
  }
  mafile->tete=NULL;
  mafile->queue=NULL;

  File2 = (file1*)malloc(sizeof(file1));
  if(!File2)
  {
    printf("erreur2\n");
    exit(-1);
  }

  File2->tete=NULL;  
  File2->queue=NULL;

  map=init_map(GridSize,GridSize,TSize);

      usedTiles=(Mybirdtab*)malloc(sizeof(Mybirdtab));
        Point=(TPoint*)malloc(sizeof(TPoint));
  flock=(birdtab*)malloc(sizeof(birdtab));

  if(!flock)
  {
    printf("Erreur \n");
    exit(-1);
  }
  flock->IDE=-1;

   steer=init_vector(0,0,0);
   vect1=init_vector(0,0,0);
   vect2=init_vector(0,0,0);
   vect5=init_vector(0,0,0);
   vect6=init_vector(0,0,0);
   vect7=init_vector(0,0,0);
   vect11=init_vector(0,0,0);

   vect3=init_vector(0,0,1);
   vect8=init_vector(0,0,1);
   vect9=init_vector(0,0,1);
   vect10=init_vector(0,0,1);
//   vect4=init_vector(0,0,1);

   averagePos=init_point(0,0,0);
   sepTmp=init_vector(0,0,0);
   cohTmp=init_point(0,0,0);
   algnTmp=init_vector(0,0,0);

  // Bird *bird;
  // for(int i=1;i<=1;i++)main
  //   for(int j=1;j<=2;j++)
  //    {

  //     bird=init_bird(init_point(i,0,j),init_vector(1,0,0),1);

  //     inserer_tab_bird(flock,bird);
  //    }

}




void draw()
{
      glTranslatef(0.0,0.0,-50.0);
      glRotatef(90,1.0,0.0,0.0);


      
     

      glColor3f(0.0, 0.0, 1.0);
      for(int i=0;i<flock->IDE+1;i++) draw_bird(flock->tab[i]);

      if(map != NULL) draw_map(map);

}









GtkWidget *fenetre,*glarea,*boiteh,*boitev,*boitev2,*bouton,*frame,*pEntry,*pLabel,*pButton,*pVBox,*label;
GSList *groupe;
int listeAttributs[] = {
  GDK_GL_RGBA,
  GDK_GL_DOUBLEBUFFER,
  GDK_GL_DEPTH_SIZE, 1,
  GDK_GL_NONE
};
int theta=-35,phi=20,xprec,yprec;
// float distance=DISTANCE_INIT;

/* Parametres d'éclairage */
GLfloat L0pos[]={ 0.0,2.0,1.0};
GLfloat L0dif[]={ 0.3,0.3,0.8};
GLfloat L1pos[]={ 2.0,2.0,2.0};
GLfloat L1dif[]={ 0.5,0.5,0.5};
GLfloat Mspec[]={0.5,0.5,0.5};
GLfloat Mshiny=50;


/* Prototypes des fonctions */
void creeInterface();
void affichage(GtkWidget *widget,GdkEventExpose *evenement);
void initGlarea(GtkWidget* widget);
void redimGlarea(GtkWidget* widget, GdkEventConfigure* evenement);
// void mouvementSouris(GtkWidget* widget, GdkEventMotion* evenement);
void rappelMap(GtkWidget* widget, GdkEventAny* evenement);
void modeAffichage(GtkWidget* widget, char *chaine);
void basculeEclairage(GtkWidget* widget,gpointer donnee);
void basculeFacesArrieres(GtkWidget *widget,gpointer donnee);
static gboolean idle(gpointer data);
gint glarea_button_press(GtkWidget *widget, GdkEventButton *event);
gint		glarea_motion_notify(GtkWidget *widget, GdkEventMotion *event);
void clear_all();


void saisir_plusieurs(GtkWidget *pButton, gpointer data);
void inserer_file_id(file1* mafile,int x);
int rechercher_tab(double x,double y);
int suppression_pos(int ind);
void ajouter_oiseau(Bird *b);
void play_idle();
void init_bird_file();
void on_activate_entry(GtkWidget *pEntry, gpointer data);
void on_copier_button(GtkWidget *pButton, gpointer data);

// srand(GTime(NULL));
void ajouter_oiseau(Bird *b)
{
    // Bird* b=init_bird(init_point(x,0.0,y), init_vector(1,0,0),1);
    inserer_tab_bird(flock,b);
    gtk_widget_queue_draw(GTK_WIDGET(glarea));
}

gint glarea_button_press(GtkWidget *widget, GdkEventButton *event)
{
	gint	return_status=TRUE;

  switch(event->type) {

    case GDK_2BUTTON_PRESS:
    if(event->button==1)
    {
        double xx=(double)event->x*40/800.0;
        double zz=(double)event->y*40/800.0;
        Bird* b=init_bird(init_point(xx,0.0,zz), init_vector(1,0,0),1);
        b->ColorIset=0;
        inserer_tab_bird(flock,b);
        gtk_widget_queue_draw(GTK_WIDGET(glarea));
    }
      break;

    case GDK_3BUTTON_PRESS:
      g_print("Mouse button %d triple-click at coordinates (%lf,%lf)\n",
              event->button,event->x,event->y);
      break;

    case GDK_BUTTON_PRESS:
    if(event->button==3) {tmp=1; g_print("%d\n",tmp); init_bird_file();}

    // if(event->button==2) {tmp2=1; }

    if(event->button==4)
    {
        double x=(double)event->x*40/800.0;
        double z=(double)event->y*40/800.0;
        Bird* b=init_bird(init_point(x,0.0,z), init_vector(1,0,0),1);
        b->ColorIset=0;
        inserer_tab_bird(flock,b);
        gtk_widget_queue_draw(GTK_WIDGET(glarea));
    }
    else if(event->button==1)
    {
      int rep=0;
      int ind = rechercher_tab((double)event->x*40/800.0,(double)event->y*40/800.0);
      if(ind!=-1)
      {
        rep= suppression_pos(ind);
        if(rep) printf("L'élement à été bien supprimé\n");
      }
          
    }

      break;

    default:
      g_print("Unknown button press event\n");
      return_status=FALSE;
      
      break;


  }

  return(return_status);
}
int rechercher_tab(double x,double y)
{
  for(i=0;i<flock->IDE+1;i++)
    if(((int)flock->tab[i]->pos->x == (int)x) && ((int)flock->tab[i]->pos->z == (int)y))
      return i;
  return -1;
}

/**
suppression_pos(Matab *Montab,int p):Cette fonction permet de supprimer un element
du tableau dans une position donne
Entree :
   int indice: ce varible me permet de circuler le tableau
Sortie: Aucun
**/
int suppression_pos(int ind)
{
    int indice;
    if((ind>(flock->IDE))||(ind<0))
       return((int)0);
    for(indice=ind;indice<(flock->IDE);indice++)
        flock->tab[indice]=flock->tab[indice+1];
     flock->IDE--;
     return((int)1);
}

void inserer_file(file* mafile,double x,double y)
{       

  lis* cellule;
  cellule=(lis*)malloc(sizeof(lis));
  if(!cellule)
  {
    printf("ERREUR 2");
    exit(-1);
  }
  cellule->x=x;
  cellule->y=y;
  cellule->svt=NULL;
  if(!mafile->tete)
  {
     mafile->tete=cellule;
     mafile->queue=cellule;
  }
  else
  {
    mafile->queue->svt=cellule;
    mafile->queue=cellule;
  }
  
 
}



gint glarea_motion_notify(GtkWidget *widget, GdkEventMotion *event)
{
	gint			x, y;
	GdkRectangle		area;
	GdkModifierType		state;
  


  if (event->is_hint) {
    gdk_window_get_pointer(event->window, &x, &y, &state);
  } else {
    x = (gint)event->x;
    y = (gint)event->y;
    state = (GdkModifierType)event->state;
  }

  if(state & GDK_BUTTON1_MASK) {

    double xx=(double)event->x*40/800.0;
    double yy=(double)event->y*40/800.0;
    /* ... GLISSEMENT... */
    inserer_file(mafile,xx,yy);
    g_print("Mouse motion button 1 at coordinates (%d,%d)\n",x,y);
    g_print("Mouse motion button 1 at coordinates (%lf,%lf)\n",xx,yy);


    /* ... Orientation has changed, redraw mesh ... */
    gtk_widget_draw(glarea, (GdkRectangle *)NULL);
  }

  if(state & GDK_BUTTON2_MASK) {
    /* ... Zooming drag ... */
    g_print("Mouse motion button 2 at coordinates (%d,%d)\n",x,y);

    /* ... Zoom has changed, redraw mesh ... */
    gtk_widget_draw(glarea, (GdkRectangle *)NULL);
  }

  if(state & GDK_BUTTON3_MASK) {
    /* ... 3rd button drag ... */
    g_print("Mouse motion button 3 at coordinates (%d,%d)\n",x,y);

    /* ... Zoom has changed, redraw mesh ... */
    gtk_widget_draw(glarea, (GdkRectangle *)NULL);
  }

  return TRUE;
}

int main(int argc,char **argv)
{
  gtk_init(&argc,&argv);
  init_drawing();
  creeInterface();


  /* affichage de l'interface */
  gtk_widget_show_all(fenetre);

  gtk_main();
  return 0;
}

void pause_idle()
{
  if(id)
  {
   
    file1* tmp=File2;
    while(tmp->tete)
    {
      printf("remove\n");
      g_source_remove(tmp->tete->ind);
      tmp->tete=tmp->tete->svt; 
      nbspeed++;
    }
    File2->queue=NULL;

    id=0;
    
  }
  
}
void play_idle()
{
  if(start) {start =0; g_source_remove(id); id=0;};
  if(!id)
  {
    play=1;
    id=g_timeout_add(5, idle, (void*)glarea);
        inserer_file_id(File2,id);
        for(int i=0;i<nbspeed-1;i++)
          inserer_file_id(File2,g_timeout_add(5, idle, (void*)glarea));
          nbspeed=0;
  }
    
  
    
}
void inserer_file_id(file1* mafile,int x)
{       

  speedFile* cellule;
  cellule=(speedFile*)malloc(sizeof(speedFile));
  if(!cellule)
  {
    printf("ERREUR 2");
    exit(-1);
  }
  cellule->ind=x;
  cellule->svt=NULL;

  if(!mafile->tete)
  {
     mafile->tete=cellule;
     mafile->queue=cellule;
  }
  else
  {
    mafile->queue->svt=cellule;
    mafile->queue=cellule;
  }
  
 
}

void speed_up()
{
  if(id)
  {
      printf("speed_up\n");
      inserer_file_id(File2,g_timeout_add(5, idle, (void*)glarea));
  }

}

void speed_down()
{
      
      if(File2->tete)
      {
        printf("speed_down\n");
        g_source_remove(File2->tete->ind);
        if(File2->tete->svt==NULL) id=0;
        File2->tete=File2->tete->svt;
        
      }
      

}


void clear_all()
{
  flock->IDE=-1;
}

/* Fonction de creation de l'interface graphique */
void creeInterface()
{
  /* la fenetre principale */
  fenetre=gtk_window_new(GTK_WINDOW_TOPLEVEL);
  gtk_window_set_title(GTK_WINDOW(fenetre),"Simulation du troupeau :    EL MARDI AYOUB    &    AL KADERI HOUSSAM-EDDINE");
  gtk_signal_connect(GTK_OBJECT(fenetre),"destroy",GTK_SIGNAL_FUNC(gtk_main_quit),NULL);
  gtk_signal_connect(GTK_OBJECT(fenetre),"delete_event",GTK_SIGNAL_FUNC(gtk_main_quit),NULL);

  
  /* un boite verticale */
  boiteh=gtk_hbox_new(FALSE,0);
  gtk_container_add(GTK_CONTAINER(fenetre),boiteh);


  /* une glarea */
  if(gdk_gl_query() == FALSE) {
    fprintf(stderr,"Impossible d'utiliser OpenGL\n");
    exit(1);
  }
  glarea = gtk_gl_area_new(listeAttributs);
  gtk_widget_set_events(GTK_WIDGET(glarea),
                        GDK_EXPOSURE_MASK|
                        GDK_BUTTON_PRESS_MASK|
			GDK_BUTTON_RELEASE_MASK|
			GDK_POINTER_MOTION_MASK|
                        GDK_POINTER_MOTION_HINT_MASK);
  gtk_widget_set_usize(GTK_WIDGET(glarea),800,800);
  gtk_box_pack_start(GTK_BOX(boiteh),glarea,TRUE,TRUE,0);
  gtk_signal_connect (GTK_OBJECT(glarea), "realize",
                      GTK_SIGNAL_FUNC(initGlarea), NULL);
  gtk_signal_connect (GTK_OBJECT(glarea), "expose_event",
                      GTK_SIGNAL_FUNC(affichage), NULL);

  gtk_signal_connect (GTK_OBJECT(glarea), "button_press_event",
                    GTK_SIGNAL_FUNC(glarea_button_press), NULL);  

gtk_signal_connect (GTK_OBJECT(glarea), "motion_notify_event",
                    GTK_SIGNAL_FUNC(glarea_motion_notify), NULL);  

  gtk_signal_connect (GTK_OBJECT(glarea), "configure_event",
                      GTK_SIGNAL_FUNC(redimGlarea), NULL);
  // gtk_signal_connect (GTK_OBJECT(glarea), "motion_notify_event  // while(mafile->tete)
  // {
  //   g_print("[%lf]---- [%lf]\n",mafile->tete->x,mafile->tete->y);
  //   mafile->tete=mafile->tete->svt;
  // }
  //                     GTK_SIGNAL_FUNC(mouvementSouris), NULL);  // while(mafile->tete)
  // {
  //   g_print("[%lf]---- [%lf]\n",mafile->tete->x,mafile->tete->y);
  //   mafile->tete=mafile->tete->svt;
  // }
  gtk_signal_connect (GTK_OBJECT(glarea), "map_event",
                      GTK_SIGNAL_FUNC(rappelMap), NULL);

  /* une boite verticale */
  boitev=gtk_vbox_new(FALSE,0);
  gtk_box_pack_start(GTK_BOX(boiteh),boitev,FALSE,FALSE,0);

  // /* une frame */
  // frame=gtk_frame_new("Mode d'affichage");
  // gtk_box_pack_start(GTK_BOX(boitev),frame,FALSE,FALSE,0);

  // /* un boite verticale pour la frame */
  // boitev2=gtk_vbox_new(FALSE,0);
  // gtk_container_add(GTK_CONTAINER(frame),boitev2);

  // /* un groupe de bouton radio */
  // bouton=gtk_radio_button_new_with_label(NULL,"Plein");
  // gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(bouton),TRUE);
  // gtk_box_pack_start(GTK_BOX(boitev2),bouton,FALSE,FALSE,0);
  // gtk_signal_connect(GTK_OBJECT(bouton),"toggled",GTK_SIGNAL_FUNC(modeAffichage),"a");
  
  // groupe=gtk_radio_button_group(GTK_RADIO_BUTTON(bouton));  
  // bouton=gtk_radio_button_new_with_label(groupe,"Fil");
  // gtk_box_pack_start(GTK_BOX(boitev2),bouton,FALSE,FALSE,0);
  // gtk_signal_connect(GTK_OBJECT(bouton),"toggled",GTK_SIGNAL_FUNC(modeAffichage),"b");
  

  // groupe=gtk_radio_button_group(GTK_RADIO_BUTTON(bouton));  
  // bouton=gtk_radio_button_new_with_label(groupe,"Points");
  // gtk_box_pack_start(GTK_BOX(boitev2),bouton,FALSE,FALSE,0);
  // gtk_signal_connect(GTK_OBJECT(bouton),"toggled",GTK_SIGNAL_FUNC(modeAffichage),"c");


  

  // /* un bouton pour bascule pour l'éclairage */
  bouton=gtk_check_button_new_with_label("Eclairage");
  gtk_box_pack_start(GTK_BOX(boitev),bouton,FALSE,FALSE,0);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(bouton),TRUE);
  gtk_signal_connect(GTK_OBJECT(bouton),"event",GTK_SIGNAL_FUNC(basculeEclairage),NULL);



  // /* un bouton pour le masquage des faces arrieres*/
  bouton=gtk_check_button_new_with_label("Masquage ");
  gtk_box_pack_start(GTK_BOX(boitev),bouton,FALSE,FALSE,0);
  gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(bouton),TRUE);
  gtk_signal_connect(GTK_OBJECT(bouton),"event",GTK_SIGNAL_FUNC(basculeFacesArrieres),NULL);




    //   glDisable(GL_CULL_FACE);
    //   glDisable(GL_LIGHTING);

    // gtk_widget_queue_draw(GTK_WIDGET(glarea));
  
  
  /* un bouton Quitter */
  bouton=gtk_button_new_with_label("Quitter");
  gtk_box_pack_end(GTK_BOX(boitev),bouton,FALSE,FALSE,0);
  gtk_signal_connect(GTK_OBJECT(bouton),"clicked",GTK_SIGNAL_FUNC(gtk_main_quit),NULL);

  frame=gtk_frame_new("Fonctionnalitees");
  gtk_box_pack_start(GTK_BOX(boitev),frame,FALSE,FALSE,0);
    /* un bouton Quitter */
 

  boitev2=gtk_vbox_new(FALSE,0);
  gtk_container_add(GTK_CONTAINER(frame),boitev2);


  bouton=gtk_button_new_with_label("Pause");
  gtk_box_pack_end(GTK_BOX(boitev2),bouton,FALSE,FALSE,0);
  gtk_signal_connect(GTK_OBJECT(bouton),"clicked",GTK_SIGNAL_FUNC(pause_idle),NULL);
  
  bouton=gtk_button_new_with_label("Play");
  gtk_box_pack_end(GTK_BOX(boitev2),bouton,FALSE,FALSE,0);
  gtk_signal_connect(GTK_OBJECT(bouton),"clicked",GTK_SIGNAL_FUNC(play_idle),NULL);

   bouton=gtk_button_new_with_label("Effacer tout");
   gtk_box_pack_end(GTK_BOX(boitev2),bouton,FALSE,FALSE,0);
   gtk_signal_connect(GTK_OBJECT(bouton),"clicked",GTK_SIGNAL_FUNC(clear_all),NULL);

  frame=gtk_frame_new("Vitesse");
  gtk_box_pack_start(GTK_BOX(boitev),frame,FALSE,FALSE,0);
    /* un bouton Quitter */
 

  boitev2=gtk_vbox_new(FALSE,0);
  gtk_container_add(GTK_CONTAINER(frame),boitev2);

   bouton=gtk_button_new_with_label("+");
   gtk_box_pack_end(GTK_BOX(boitev2),bouton,FALSE,FALSE,0);
   gtk_signal_connect(GTK_OBJECT(bouton),"clicked",GTK_SIGNAL_FUNC(speed_up),NULL);

   bouton=gtk_button_new_with_label("-");
   gtk_box_pack_end(GTK_BOX(boitev2),bouton,FALSE,FALSE,0);
   gtk_signal_connect(GTK_OBJECT(bouton),"clicked",GTK_SIGNAL_FUNC(speed_down),NULL);



     /* une frame */
    frame=gtk_frame_new("Saisie");
    gtk_box_pack_start(GTK_BOX(boitev),frame,FALSE,FALSE,0);

    pVBox = gtk_vbox_new(TRUE, 0);
    gtk_container_add(GTK_CONTAINER(frame), pVBox);

    pEntry = gtk_entry_new();
    label = gtk_label_new ("Position X (0-40)");

    // gtk_entry_set_placeholder_text(GTK_ENTRY (text), "x : 0 - 40");
    /* Insertion du GtkEntry dans la GtkVBox */
    gtk_box_pack_start(GTK_BOX(pVBox), label, TRUE, FALSE, 0);
    gtk_box_pack_start(GTK_BOX(pVBox), pEntry, TRUE, FALSE, 0);
    
    pEntry = gtk_entry_new();
    label = gtk_label_new ("Position Y (0-40)");
    gtk_box_pack_start(GTK_BOX(pVBox), label, TRUE, FALSE, 0);
    /* Insertion du GtkEntry dans la GtkVBox */
    gtk_box_pack_start(GTK_BOX(pVBox), pEntry, TRUE, FALSE, 0);

    pEntry = gtk_entry_new();
    label = gtk_label_new ("Age  ");
    gtk_box_pack_start(GTK_BOX(pVBox), label, TRUE, FALSE, 0);
    /* Insertion du GtkEntry dans la GtkVBox */
    gtk_box_pack_start(GTK_BOX(pVBox), pEntry, TRUE, FALSE, 0);

    
    pEntry = gtk_entry_new();
    label = gtk_label_new ("Rouge (0.0 - 1.0)");
    gtk_box_pack_start(GTK_BOX(pVBox), label, TRUE, FALSE, 0);
    /* Insertion du GtkEntry dans la GtkVBox */
    gtk_box_pack_start(GTK_BOX(pVBox), pEntry, TRUE, FALSE, 0);

    pEntry = gtk_entry_new();
    label = gtk_label_new ("Vert (0.0 - 1.0)");
    gtk_box_pack_start(GTK_BOX(pVBox), label, TRUE, FALSE, 0);
    /* Insertion du GtkEntry dans la GtkVBox */
    gtk_box_pack_start(GTK_BOX(pVBox), pEntry, TRUE, FALSE, 0);

    pEntry = gtk_entry_new();
    label = gtk_label_new ("Bleu (0.0 - 1.0)");
    gtk_box_pack_start(GTK_BOX(pVBox), label, TRUE, FALSE, 0);
    /* Insertion du GtkEntry dans la GtkVBox */
    gtk_box_pack_start(GTK_BOX(pVBox), pEntry, TRUE, FALSE, 0);





    // pLabel = gtk_label_new(NULL);
    // gtk_box_pack_start(GTK_BOX(pVBox), pLabel, TRUE, FALSE, 0);

    pButton = gtk_button_new_with_label("Valider");
    gtk_box_pack_start(GTK_BOX(pVBox), pButton, TRUE, FALSE, 0);

    // gtk_signal_connect(GTK_OBJECT(pEntry), "activate", GTK_SIGNAL_FUNC(on_activate_entry), (GtkWidget*) pLabel);


    gtk_signal_connect(GTK_OBJECT(pButton), "clicked", GTK_SIGNAL_FUNC(on_copier_button), (GtkWidget*) pVBox);

  // gtk_signal_connect(GTK_OBJECT(bouton),"clicked",GTK_SIGNAL_FUNC(gtk_main_quit),NULL);



       /* une frame */
    frame=gtk_frame_new("Saisir plusieurs");
    gtk_box_pack_start(GTK_BOX(boitev),frame,FALSE,FALSE,0);

    pVBox = gtk_vbox_new(TRUE, 0);
    gtk_container_add(GTK_CONTAINER(frame), pVBox);

    pEntry = gtk_entry_new();
    /* Insertion du GtkEntry dans la GtkVBox */
    gtk_box_pack_start(GTK_BOX(pVBox), pEntry, TRUE, FALSE, 0);

    pButton = gtk_button_new_with_label("OK");
    gtk_box_pack_start(GTK_BOX(pVBox), pButton, TRUE, FALSE, 0);

    // gtk_signal_connect(GTK_OBJECT(pEntry), "activate", GTK_SIGNAL_FUNC(on_activate_entry), (GtkWidget*) pLabel);


    
    gtk_signal_connect(GTK_OBJECT(pButton), "clicked", GTK_SIGNAL_FUNC(saisir_plusieurs), (GtkWidget*) pVBox);

}


void on_activate_entry(GtkWidget *pEntry, gpointer data)
{
    const gchar *sText;
 
    /* Recuperation du texte contenu dans le GtkEntry */
    sText = gtk_entry_get_text(GTK_ENTRY(pEntry));
    
    /* Modification du texte contenu dans le GtkLabel */
    gtk_label_set_text(GTK_LABEL((GtkWidget*)data), sText);
}

void saisir_plusieurs(GtkWidget *pButton, gpointer data)
{
  GtkWidget *pTempEntry;
    // GtkWidget *pTempLabel;
    GList *pList;
    const gchar *sText;
 
    /* Recuperation de la liste des elements que contient la GtkVBox */
    pList = gtk_container_children(GTK_CONTAINER((GtkWidget*)data));
    
 
    /* Le premier element est le GtkEntry */
    pTempEntry = GTK_WIDGET(pList->data);
    sText = gtk_entry_get_text(GTK_ENTRY(pTempEntry));
    int nbr= atoi(sText);
    // sText = gtk_entry_get_text(GTK_ENTRY(pTempEntry2));
    // int y = atoi(sText);
    for(int i=0;i<nbr;i++)
    {
      Bird* b=init_bird(init_point(rand()%40,0.0,rand()%40), init_vector(1,0,0),1);
      b->ColorIset = 0;
      ajouter_oiseau(b);
    }
       

     g_list_free(pList);
}

void on_copier_button(GtkWidget *pButton, gpointer data)
{
    GtkWidget *pTempEntry,*pTempEntry2,*pTempEntry3,*pTempEntry4,*pTempEntry5,*pTempEntry6;
    GtkWidget *pTempLabel;
    GList *pList;
    const gchar *sText;
 
    /* Recuperation de la liste des elements que contient la GtkVBox */
    pList = gtk_container_children(GTK_CONTAINER((GtkWidget*)data));
    
    
        /* Le premier element est le GtkEntry */
      pList = g_list_next(pList);
      pTempEntry = GTK_WIDGET(pList->data);
      pList = g_list_next(pList);
      pList = g_list_next(pList);
      pTempEntry2 = GTK_WIDGET(pList->data);
      pList = g_list_next(pList);
      pList = g_list_next(pList);
      pTempEntry3 = GTK_WIDGET(pList->data);
      pList = g_list_next(pList);
      pList = g_list_next(pList);
      pTempEntry4 = GTK_WIDGET(pList->data);
      pList = g_list_next(pList);
      pList = g_list_next(pList);
      pTempEntry5 = GTK_WIDGET(pList->data);
      pList = g_list_next(pList);
      pList = g_list_next(pList);
      pTempEntry6 = GTK_WIDGET(pList->data);
    

 
    /* Passage a l element suivant : le GtkButton */

    
 
    /* Passage a l element suivant : le GtkLabel */

    // 
 
    /* Cet element est le GtkLabel */

    // pTempLabel = GTK_WIDGET(pList->data);
// printf("%p\n",pTempEntry);
          /* Recuperation du texte contenu dans le GtkEntry */
    sText = gtk_entry_get_text(GTK_ENTRY(pTempEntry));
    int x= atoi(sText);
   
    sText = gtk_entry_get_text(GTK_ENTRY(pTempEntry2));
    int y = atoi(sText);

    Bird* b=init_bird(init_point(x,0.0,y), init_vector(1,0,0),1);

    sText = gtk_entry_get_text(GTK_ENTRY(pTempEntry3));
    b->age= atof(sText);

    sText = gtk_entry_get_text(GTK_ENTRY(pTempEntry4));
    b->red = atof(sText);

    sText = gtk_entry_get_text(GTK_ENTRY(pTempEntry5));
    b->green= atof(sText);

    sText = gtk_entry_get_text(GTK_ENTRY(pTempEntry6));
    b->blue= atof(sText);

    b->ColorIset =1;
    
    ajouter_oiseau(b);
    /* Modification du texte contenu dans le GtkLabel */

    // gtk_label_set_text(GTK_LABEL(pTempLabel), sText);
 
    /* Liberation de la memoire utilisee par la liste */
    g_list_free(pList);
    
    

}
Bird* bird;
/* Fonction de rappel pour l'affichage */
void affichage(GtkWidget *widget,GdkEventExpose *evenement)
{
  if (evenement->count > 0)
    return;

  if (gtk_gl_area_make_current(GTK_GL_AREA(widget))) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();


    //     gtk_gl_area_make_current(GTK_GL_AREA(glarea));
    // glEnable(GL_CULL_FACE);
    // gtk_widget_queue_draw(GTK_WIDGET(glarea));
    // glEnable(GL_LIGHTING);
    // gtk_widget_queue_draw(GTK_WIDGET(glarea));
    draw();





    gtk_gl_area_swapbuffers(GTK_GL_AREA(widget));
  }
}


/* fonction d'initialisation d'OpenGL*/
void initGlarea(GtkWidget* widget)
{
  if (gtk_gl_area_make_current(GTK_GL_AREA(widget))) {
    // glClearColor(0.4,0.4,0.4,1.0);
    glClearColor(1.0,1.0,1.0,0.0);
    glColor3d(1,0,0);
    glPointSize(4.0);
    glLineWidth(2.0);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);


  int l = widget->allocation.width;
  int h = widget->allocation.height;

  glViewport(0,-65,l, h+130);

    /* Mise en place de la perspective */
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();  
    // gluPerspective(45.0,1.0,0.1,DISTANCE_MAX+2);
    gluPerspective(50.0, (float) l / (float) (h+130), 1.0, 100.0);

    glMatrixMode(GL_MODELVIEW);


    /* Parametrage de l'éclairage */
    glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER,GL_TRUE);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1);
    glLightfv(GL_LIGHT0,GL_DIFFUSE,L0dif);
    glLightfv(GL_LIGHT0,GL_SPECULAR,L0dif);
    glLightfv(GL_LIGHT1,GL_DIFFUSE,L1dif);
    glLightfv(GL_LIGHT1,GL_SPECULAR,L1dif); 
    
    /* Paramétrage du matériau */
    glMaterialfv(GL_FRONT_AND_BACK,GL_SPECULAR,Mspec);
    glMaterialf(GL_FRONT_AND_BACK,GL_SHININESS,Mshiny);




    id=g_timeout_add(5, idle, (void*)widget);

    start=1;
  }
}
// void tracer_chemin(Bird *b)
// {
//   while(mafile->tete->svt)
//   {
//     double x = (double)(((mafile->tete->x)*40)/800.0);
//     double y = (double)(((mafile->tete->y)*40)/800.0);
    
//     mafile->tete=mafile->tete->svt;
//   }
//   tmp=0;
// }


void init_bird_file()
{ 
  if(mafile->tete)
  if(mafile->tete->svt)
  {
      
      double x = mafile->tete->x;
      double y = mafile->tete->y;
      // g_print("[%lf]---- [%lf]\n",x,y);
      bird=init_bird(init_point(x,0.0,y), init_vector(1,0,0),1);
      bird->ColorIset=0;
      if(mafile->tete->svt->x > x)
      {
        if(mafile->tete->svt->y > y)
          bird->angle=45;
        else if(mafile->tete->svt->y == y)
          bird->angle=90;
        else if(mafile->tete->svt->y < y)
          bird->angle=135;
      }
      else if(mafile->tete->svt->x == x)
      {
          if(mafile->tete->svt->y > y)
          bird->angle=0;
          
        else if(mafile->tete->svt->y < y)
          bird->angle=180;
      }
      else if(mafile->tete->svt->x < x)
      {
        if(mafile->tete->svt->y > y)
          bird->angle=-45;
        else if(mafile->tete->svt->y == y)
          bird->angle=-90;
        else if(mafile->tete->svt->y < y)
          bird->angle=-135;
      }
      
        

      glTranslatef(0.0,0.0,-50.0);
      glRotatef(90,1.0,0.0,0.0);
      glColor3f(0.0, 0.0, 1.0);
      draw_bird(bird);
      inserer_tab_bird(flock,bird);
      gtk_widget_queue_draw(GTK_WIDGET(glarea));
      mafile->tete=mafile->tete->svt;
  }

}

static gboolean idle(gpointer data) {
if(play)
    flockMovement();
  if(tmp)
  {
      if(mafile->tete)
      if(mafile->tete->svt)
      {
          
          bird->pos->x = mafile->tete->x;
          bird->pos->z = mafile->tete->y;
          if(mafile->tete->svt->x > bird->pos->x)
          {
          if(mafile->tete->svt->y > bird->pos->z)
            bird->angle=45;
          else if(mafile->tete->svt->y == bird->pos->z)
            bird->angle=90;
          else if(mafile->tete->svt->y < bird->pos->z)
            bird->angle=135;
          }
          else if(mafile->tete->svt->x == bird->pos->x)
          {
            if(mafile->tete->svt->y > bird->pos->z)
            bird->angle=0;
            
          else if(mafile->tete->svt->y < bird->pos->z)
            bird->angle=180;
          }
          else if(mafile->tete->svt->x < bird->pos->x)
          { 
          if(mafile->tete->svt->y > bird->pos->z)
            bird->angle=-45;
          else if(mafile->tete->svt->y == bird->pos->z)
            bird->angle=-90;
          else if(mafile->tete->svt->y < bird->pos->z)
            bird->angle=-135;
          }
        glTranslatef(0.0,0.0,-50.0);
        glRotatef(90,1.0,0.0,0.0);
        glColor3f(0.0, 0.0, 1.0);
        draw_bird(bird);
        gtk_widget_queue_draw(GTK_WIDGET(glarea));
        mafile->tete=mafile->tete->svt;
    }
    else tmp=0;

  }
    
	gtk_widget_queue_draw(GTK_WIDGET(data));
	return TRUE;

}






/* Fonction de rappel pour le redimensionnement du widget OpenGL  */
void redimGlarea(GtkWidget* widget, GdkEventConfigure* evenement) 
{
  int l = widget->allocation.width;
  int h = widget->allocation.height;
  if (gtk_gl_area_make_current(GTK_GL_AREA(widget))) {
    // if (l<h)
    //   glViewport(0,(h-l)/2,l,l);
    // else 
    //   glViewport((l-h)/2,0,h,h);

  glViewport(0,0,l, h-10);
  //make changes to the projection matrix
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluPerspective(50.0, (float) l / (float) (h+80), 1.0, 100.0);


  }
}



/* fonction de rappel pour les mouvements de la souris */
// void mouvementSouris(GtkWidget* widget, GdkEventMotion* evenement)
// {
//   int x,y;
//   GdkModifierType etat;
//   if (evenement->is_hint)
//     gdk_window_get_pointer(evenement->window,&x,&y,&etat);
//   else {
//     x=evenement->x;
//     y=evenement->y;
//     etat=evenement->state;
//   }

//   /* bouton gauche = rotation */
//   if (etat & GDK_BUTTON1_MASK){
//     theta+=x-xprec;
//     phi+=y-yprec;
//     gtk_widget_queue_draw(widget);
//   }

//   /* bouton droit = Zoom */
//   if (etat & GDK_BUTTON3_MASK){
//     distance+=((float)(y-yprec))/20.0;
//     if (distance<DISTANCE_MIN)
//       distance=DISTANCE_MIN;
//     if (distance>DISTANCE_MAX)
//       distance=DISTANCE_MAX;
//     gtk_widget_queue_draw(widget);
//   }
//   xprec=x;yprec=y;
// }



/* Fonction de rappel pour le mapping du widget */
void rappelMap(GtkWidget* widget, GdkEventAny* evenement)
{
  int x,y;
  gdk_window_get_pointer(glarea->window,&x,&y,NULL);
  xprec=x;
  yprec=y;
}

/* Fonction de bascule du mode d'affichage */
void modeAffichage(GtkWidget* widget, char *chaine)
{
  if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget)))
    if (gtk_gl_area_make_current(GTK_GL_AREA(glarea))) {
      switch (chaine[0]) {
      case 'a':
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
	break;
      case 'b':
	glPolygonMode(GL_FRONT_AND_BACK,GL_LINE);
	break;
      case 'c':
	glPolygonMode(GL_FRONT_AND_BACK,GL_POINT);
	break;
      }
      gtk_widget_queue_draw(GTK_WIDGET(glarea));
    }
}



/* Fonction de bascule du mode d'eclairage */
void basculeEclairage(GtkWidget* widget,gpointer donnee)
{
  // if (gtk_gl_area_make_current(GTK_GL_AREA(glarea))) {
    // if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget)))
      // glEnable(GL_LIGHTING);
    // else
      glDisable(GL_LIGHTING);
    gtk_widget_queue_draw(GTK_WIDGET(glarea));
  // }
}



/* Fonction de bacule pour le masquage des faces arrieres */
void basculeFacesArrieres(GtkWidget *widget,gpointer donnee)
{
  // if (gtk_gl_area_make_current(GTK_GL_AREA(glarea))) {
    // if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(widget)))
      // glEnable(GL_CULL_FACE);
    // else
      glDisable(GL_CULL_FACE);
    gtk_widget_queue_draw(GTK_WIDGET(glarea));
  // }
}