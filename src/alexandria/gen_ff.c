#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "typedefs.h"
#include "atomprop.h"
#include "poldata.h"
#include "poldata_xml.h"

int main(int argc,char *argv[]) 
{
  gmx_poldata_t pd;
  char *smname,*elem,*desc,*gt_type,*gt_old,*charge;
  char *neighbors,*vdwparams,*geometry;
  int numbonds,atomnumber;
  int btp=1,atp=1,dtp=1;
  char *ai,*aj,*ak,*al,*params;
  double polarizability,sig_pol;
  real mass;
  gmx_atomprop_t aps;
  FILE *fp;
  
  aps = gmx_atomprop_init();
  pd  = gmx_poldata_read(NULL,aps);
  fp = fopen("forcefield.itp","w");
  fprintf(fp,"; This file is generated from gentop.dat by %s\n",argv[0]);
  fprintf(fp,"; Do not edit this file!\n");
  fprintf(fp,"; This is the force field file for the Alexandria FF\n");
  fprintf(fp,"; Paul J. van Maaren and David van der Spoel\n\n");
  fprintf(fp,"[ defaults ]\n");
  fprintf(fp,"; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n");
  fprintf(fp,"1               1               yes             1       1\n\n");

  fprintf(fp,"[ atomtypes ]\n");
  fprintf(fp,"; name      at.num  mass     charge ptype  repulsion  dispersion\n");

  gt_old = NULL;
  while (1 == gmx_poldata_get_atype(pd,&elem,&desc,&gt_type,
                                    NULL,&charge,&polarizability,&sig_pol,
                                    &vdwparams))
    {
      if (gmx_atomprop_query(aps,epropMass,"",elem,&mass)) 
	{
	  atomnumber = gmx_atomprop_atomnumber(aps,elem);
	  if ((NULL == gt_old) || (strcmp(gt_old,gt_type) != 0)) {
	    fprintf(fp,"%5s   %3d  %12.6f  %10.4f  A  %-s\n",
		    gt_type,atomnumber,mass,0.0,vdwparams);
	    fprintf(fp,"%5ss  %3d  %12.6f  %10.4f  S  0  0\n",
		    gt_type,0,0.0,0.0);
	  }
	  gt_old = gt_type;
	}
    }
  /* Bondtypes */
  fprintf(fp,"\n[ bondtypes ]\n");
  fprintf(fp,"; ; i    j  func       parameters\n");
  while (1 == gmx_poldata_get_bond(pd,&ai,&aj,NULL,NULL,NULL,&params)) {
    fprintf(fp,"%-5s  %-5s   %d  %s\n",ai,aj,btp,params);
  }

  /* Angletypes */
  fprintf(fp,"\n[ angletypes ]\n");
  fprintf(fp,"; ; i    j   k  func       parameters\n");
  while (1 == gmx_poldata_get_angle(pd,&ai,&aj,&ak,NULL,NULL,&params)) {
    fprintf(fp,"%-5s  %-5s  %-5s  %d  %s\n",ai,aj,ak,atp,params);
  }
  
  /* Dihedraltypes */
  fprintf(fp,"\n[ dihedraltypes ]\n");
  fprintf(fp,"; ; i    j   k    l  func       parameters\n");
  while (1 == gmx_poldata_get_dihedral(pd,&ai,&aj,&ak,&al,NULL,NULL,&params)) {
    fprintf(fp,"%-5s  %-5s  %-5s  %-5s  %d  %s\n",ai,aj,ak,al,dtp,params);
  }
  
  fclose(fp);
  
  return 0;
}