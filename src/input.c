#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include<stdlib.h>
#include "input.h"
#include "types.h"
#include "utilities.h"

/* helper function: read a line and then return
   the first string with whitespace stripped off */
int get_a_line(FILE *fp, char *buf)
{
    char tmp[BLEN], *ptr;

    /* read a line and cut of comments and blanks */
    if (fgets(tmp,BLEN,fp)) {
        int i;

        ptr=strchr(tmp,'#');
        if (ptr) *ptr= '\0';
        i=strlen(tmp); --i;
        while(isspace(tmp[i])) {
            tmp[i]='\0';
            --i;
        }
        ptr=tmp;
        while(isspace(*ptr)) {++ptr;}
        i=strlen(ptr);
        strcpy(buf,tmp);
        return 0;
    } else {
        perror("problem reading input");
        return -1;
    }
    return 0;
}



void read_from_file(mdsys_t * sys, int * nprint,char * restfile, char* trajfile, char * ergfile,char * line   ){
     /* read input file */
    if(get_a_line(stdin,line)) exit(1);
    (*sys).natoms=atoi(line);
    if(get_a_line(stdin,line)) exit(1);
    (*sys).mass=atof(line);
    if(get_a_line(stdin,line)) exit(1);
    (*sys).epsilon=atof(line);
    if(get_a_line(stdin,line)) exit(1);
    (*sys).sigma=atof(line);
    if(get_a_line(stdin,line)) exit(1);
    (*sys).rcut=atof(line);
    if(get_a_line(stdin,line))  exit(1);
    (*sys).box=atof(line);
    if(get_a_line(stdin,restfile)) exit(1);
    if(get_a_line(stdin,trajfile)) exit(1);
    if(get_a_line(stdin,ergfile)) exit(1);
    if(get_a_line(stdin,line)) exit(1);
    (*sys).nsteps=atoi(line);
    if(get_a_line(stdin,line)) exit(1);
    (*sys).dt=atof(line);
    if(get_a_line(stdin,line)) exit(1);
    *nprint=atoi(line);
}





    /* read restart */
void read_restfile(char *restfile, mdsys_t *sys){
 int i;
 FILE *fp;

  fp=fopen(restfile,"r");
    if(fp) {
        for (i=0; i<sys->natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys->rx+i, sys->ry+i, sys->rz+i);
        }
        for (i=0; i<sys->natoms; ++i) {
            fscanf(fp,"%lf%lf%lf",sys->vx+i, sys->vy+i, sys->vz+i);
        }
        fclose(fp);
        azzero(sys->fx, sys->natoms);
        azzero(sys->fy, sys->natoms);
        azzero(sys->fz, sys->natoms);
    } else {
        perror("cannot read restart file");
        exit(3);
    }
}
