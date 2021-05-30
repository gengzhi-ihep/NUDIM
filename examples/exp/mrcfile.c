 #include "mrcfile.h"


 /*******************************************************************************************/
 int mrc_read_head (FILE *fin,  MrcHeader *head)
 {
   if(ftell(fin)!=0)rewind(fin);
   fread(head,1024,1,fin);

   if(!(head->cmap[0]=='M'&&head->cmap[1]=='A'&&head->cmap[2]=='P'))
   {
      printf("Error with function 'mrc_read_head()'! Warning: Not MRC format! \n");
      return -1;
   }

   return 0;
 }

 /*******************************************************************************************/
 int mrc_init_head(MrcHeader *head)
 {
   head->nx=0;
   head->ny=0;
   head->nz=0;

   head->mode=MRC_MODE_FLOAT;

   head->nxstart=0;
   head->nystart=0;
   head->nzstart=0;

   head->mx=1;
   head->my=1;
   head->mz=1;

   head->xlen=1;
   head->ylen=1;
   head->zlen=1;

   head->alpha=90;
   head->beta=90;
   head->gamma=90;

   head->mapc=1;
   head->mapr=2;
   head->maps=3;

   head->amin=0;
   head->amax=255;
   head->amean=128;

   head->ispg=0;
   head->nsymbt=0;

   head->next=0;

   head->creatid=1000;
   head->cmap[0]='M';
   head->cmap[1]='A';
   head->cmap[2]='P';

   head->stamp[0]='D';
 }

 /*******************************************************************************************/
 int mrc_read_slice(FILE *fin, MrcHeader  *head, int slcN,char axis,float *slcdata)
 {
   int psize;
   short buf_short;
   unsigned char buf_byte;
   float buf_float;
   int i,j,k;

   switch(head->mode)
   {
     case MRC_MODE_BYTE :
     psize=sizeof(unsigned char);

     break;

     case MRC_MODE_SHORT :
     psize=sizeof(short);

     break;

     case MRC_MODE_FLOAT :
     psize=sizeof(float);

     break;
   }

   switch(axis)
   {

 /***********************************X************************************/
     case 'x':
     case 'X':

       fseek(fin,1024+slcN*psize,SEEK_SET);


       switch(head->mode)
       {
         case MRC_MODE_BYTE:
           for(i=0;i<head->ny*head->nz;i++)
               {
                 fread(&buf_byte,psize,1,fin);
                 slcdata[i]=(float)buf_byte;
                 fseek(fin,(head->nx-1)*psize,SEEK_CUR);
               }

           break;

         case MRC_MODE_SHORT:
           for(i=0;i<head->ny*head->nz;i++)
              {
                 fread(&buf_short,psize,1,fin);
                 slcdata[i]=(float)(buf_short);
                 fseek(fin,(head->nx-1)*psize,SEEK_CUR);
               }

           break;

         case MRC_MODE_FLOAT:
           for(i=0;i<head->ny*head->nz;i++)
             {
             fread(&buf_float,psize,1,fin);
             slcdata[i]=buf_float;
             fseek(fin,(head->nx-1)*psize,SEEK_CUR);
             }
           break;

        }

     break;

 /***********************************Y************************************/
     case 'y':
     case 'Y':

       for(k=0;k<head->nz;k++)
       {
         fseek(fin,1024+(k*head->nx*head->ny+head->nx*slcN)*psize,SEEK_SET);


       switch(head->mode)
       {
         case MRC_MODE_BYTE:
         for(i=0;i<head->nx;i++)
               {
                 fread(&buf_byte,psize,1,fin);
                 slcdata[k*head->nx+i]=(float)buf_byte;
               }

           break;

         case MRC_MODE_SHORT:
         for(i=0;i<head->nx;i++)
              {
                 fread(&buf_short,psize,1,fin);
                 slcdata[k*head->nx+i]=(float)(buf_short);
               }

           break;

         case MRC_MODE_FLOAT:
         for(i=0;i<head->nx;i++)
             {
             fread(&buf_float,psize,1,fin);
             slcdata[k*head->nx+i]=buf_float;
             }

           break;

       }


       }
     break;

 /***********************************Z************************************/
     case 'z':
     case 'Z':
       fseek(fin,1024+slcN*head->nx*head->ny*psize,SEEK_SET);

       if(head->mode==MRC_MODE_FLOAT){
		   fread(slcdata,psize*head->nx*head->ny,1,fin);
	   }

       if(head->mode==MRC_MODE_BYTE)
       {
         for(i=0;i<head->nx*head->ny;i++)
         {
           fread(&buf_byte,psize,1,fin);
           slcdata[i]=(float)buf_byte;
         }
        }

       if(head->mode==MRC_MODE_SHORT)
       {
         for(i=0;i<head->nx*head->ny;i++)
         {
           fread(&buf_short,psize,1,fin);
           slcdata[i]=(float)buf_short;
         }
        }
     break;

     }
 }




