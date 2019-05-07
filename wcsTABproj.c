/*============================================================================
*=============================================================================
*
* wcsTABproj reorganises TAB projected data into a 2-D/3-D FITS image.
* Refer to the usage notes below.
*
*---------------------------------------------------------------------------*/

char usage[] =
"Usage: wcsTABproj [OPTION]... <infile> <outfile>\n"
"\n"
"wcsTABproj reorganises spatially TAB projected data into a 2-D/3-D FITS"
" image.\n"
"Options:\n"
"  -a<alt>         Specify an alternate coordinate representation (ignored if\n"
"                  there is only one).  Can also be specified as a 0-relative\n"
"                  index in the range 0 to 26, where alternates are sequenced\n"
"                  alphabetically following the primary representation.\n"
"  -h<header_name> Specify the name of the text WCS header to be used for the \n"
"                  output image.\n"
"  -r<resolution>  Specify the resolution of the output image (deg/px).\n"
"  -n1<pixelnumber_x> Specify the x dimension of the output image.\n"
"  -n2<pixelnumber_y> Specify the y dimension of the output image.\n";


#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include <fitsio.h>

#include <wcs.h>
#include <wcshdr.h>
#include <wcsfix.h>
#include <wcsprintf.h>
#include <getwcstab.h>

int main(int argc, char **argv)

{

  FILE  *hptr;
  int   i, j, z, alts[27], *stat = 0x0, *statout = 0x0, status, nkeyrec,
        hdutype, relax, ctrl, nreject, nwcs, nelem, *st = 0x0,
        ialt, nfound, bitpix, nwcsout, nkeys, n, nn;
  long  firstpix, ipix, lastpix, naxis, *naxes = 0x0, nside = 0, repeat,
        offset, n1 = 0, n2 = 0, n3, *naxesout = 0x0;
  double *imgcrd = 0x0, phi, *pixcrd = 0x0, theta, *world = 0x0, *worldin = 0x0,
         *imgcrdout = 0x0, *pixcrdout = 0x0, res = 0.0;

  double *pixul = 0x0, *pixur = 0x0, *pixbl = 0x0, *pixbr = 0x0, *wul= 0x0,
         *wur = 0x0, *wbl = 0x0, *wbr = 0x0, *imcrd = 0x0, lnmin, lnmax, lgmin,
         lgmax;

  char  *infile;	/* Input file.                                */
  char  *outfile;	/* Output file.                               */
  char  *alt = 0x0, *header, headfile[10], outhead[2048], headline[81],
        wcsname[72], wcsnameout[72];

  size_t nb;

  struct wcsprm *wcs, *wcsi = 0x0, *wcsout;
  fitsfile *fptr, *fptrout;

  /* Initialise control variables. */
  ctrl   = -3;
  relax  = WCSHDR_all;


  /* Parse options. */
  for (i = 1; i < argc && argv[i][0] == '-'; i++) {
    if (!argv[i][1]) break;
    switch (argv[i][1]) {
    case 'a':
      /* Select an alternate WCS. */
      alt = argv[i]+2;
      break;

    case 'h':
      /* Output WCS header. */
      sscanf(argv[i]+2,"%s",&headfile);
      break;

    case 'r':
      /* Output resolution. */
      res = (double)(atof(argv[i]+2));
      break;

    case 'n':
      /* Output dimension. */
      switch (argv[i][2]) {

      case '1':
        n1 = atof(argv[i]+3);
        break;

      case '2':
        n2 = atof(argv[i]+3);
        break;

      default:
        wcsfprintf(stderr, "%s", usage);
        return 1;
      }

      break;

    default:
      wcsfprintf(stderr, "%s", usage);
      return 1;
    }
  }
  if (i < argc) {
    infile = argv[i++];

    if (i < argc) {
      outfile = argv[i++];

      if (i < argc) {
        fprintf(stderr, "%s", usage);
        return 1;
      }
    } else {
        fprintf(stderr, "%s", usage);
        return 1;
    }
  } else {
      fprintf(stderr, "%s", usage);
      return 1;
  }
 

  /* Open the FITS file. */
  status = 0;
  if (fits_open_file(&fptr, infile, READONLY, &status)) goto fitserr;

  /* Read in the FITS header, excluding COMMENT and HISTORY keyrecords. */
  if (fits_hdr2str(fptr, 1, NULL, 0, &header, &nkeyrec, &status)) {
    goto fitserr;
  }
  if (fits_get_hdu_type(fptr, &hdutype, &status)) goto fitserr;

  /* Get the image pixel value type */
  if (fits_get_img_type(fptr, &bitpix, &status)) goto fitserr;

  /* Get the image size. */
  if (fits_read_key_lng(fptr, "NAXIS", &naxis, 0x0, &status)) goto fitserr;
  naxes = malloc(naxis*sizeof(long));
  if (fits_read_keys_lng(fptr, "NAXIS", 1, (int)naxis, naxes, &nfound,
                         &status)) goto fitserr;

  /* Interpret the WCS keywords. */
  if (hdutype == IMAGE_HDU) {
    if ((status = wcspih(header, nkeyrec, relax, ctrl, &nreject,
                    &nwcs, &wcs))) {
      wcsfprintf(stderr, "wcspih ERROR %d: %s.\n", status,
                 wcshdr_errmsg[status]);
      return 1;
    }
  } else {
    wcsfprintf(stderr, "wcsware: Invalid FITS extension type.\n");
    return 1;
  }

  free(header);

  if (nreject) {
    if (ctrl <= 3) {
      wcsprintf("\n%d WCS keyrecords were rejected.\n", nreject);
    }
  } else if (nwcs == 0) {
    if (2 < ctrl) wcsprintf("\n");
    wcsprintf("No world coordinate systems found.\n");
    fits_close_file(fptr, &status);
    return 0;

  } else if (2 < ctrl) {
    wcsprintf("\nNo invalid WCS keyrecords were found.\n");
  }

  /* Sort out alternates. */
  i = 0;
  if (alt) {

    if ('0' <= *alt && *alt <= '9') {
      if ((i = atoi(alt)) > nwcs-1) {
        wcsfprintf(stderr, "WARNING, no alternate coordinate "
          "representation \"%s\".\n", alt);
        return 1;
      }

    } else {
      wcsidx(nwcs, &wcs, alts);

      ialt = toupper(*alt);
      if (strlen(alt) > 1) {
        wcsfprintf(stderr, "WARNING, alternate specifier \"%s\" is "
          "invalid.\n", alt);
        return 1;

      } else if (*alt == ' ') {
        if (alts[0] == -1) {
          wcsfprintf(stderr, "WARNING, no primary coordinate "
            "representation.\n");
          return 1;
        }

      } else if (ialt < 'A' || ialt > 'Z') {
        wcsfprintf(stderr, "WARNING, alternate specifier \"%s\" is "
          "invalid.\n", alt);
        return 1;

      } else {
        if ((i = alts[ialt - 'A' + 1]) == -1) {
          wcsfprintf(stderr, "WARNING, no alternate coordinate "
            "representation \"%s\".\n", alt);
          return 1;
        }
      }
    }
  }
  wcsi = wcs + i;

  /* Initialize and possibly print the structs. */
  wcserr_enable(1);
  wcsprintf_set(stdout);


  /* Read -TAB arrays from the binary table extension (if necessary). */
  if (fits_read_wcstab(fptr, wcs[i].nwtb, (wtbarr *)wcs[i].wtb,
                         &status)) {
      goto fitserr;
  }

  /* Get WCSNAME out of the wcsprm struct. */
  strcpy(wcsname, wcs[i].wcsname);
  printf("%s\n",wcsname);

  /* Get the output WCS structure from ascii header */
  if (headfile) {
    if (hptr = fopen(headfile, "r")) {
      nkeys = 0;
      while (fgets(headline, 81, hptr)) {
          if ((nb = strlen(headline)) < 81) {
              for (nn=nb-1; nn<79; nn++) {
                  strncpy(headline+nn, " ", 1);
              }
              strncpy(headline+nn, "\n\0", 2);              
          }
          offset = nkeys * 80;
          strncpy(outhead+offset,headline, 80);
          nkeys++;
      }
      fclose(hptr);
      if ((status = wcspih(outhead, nkeys, relax, ctrl, &nreject,
                    &nwcsout, &wcsout))) {
          wcsfprintf(stderr, "wcspih ERROR %d: %s.\n", status,
                     wcshdr_errmsg[status]);
          return 1;
      }
    } else {
      fprintf(stderr, "ERROR: cannot open header file \"%s\".\n", headfile);
      return 1;
    }
  }


  /* Get the output WCS structure from cli */
  if (res > 0.0) {
    wcsout->naxis = wcsi->naxis;
    for (n=0; n<wcsout->naxis; n++) {
      wcsout->cdelt[n] = res;
    }
  }

  /* Get the image corners in world coordinates */
  imcrd = realloc(imcrd, naxis * sizeof(double));
  st = realloc(st, naxis * sizeof(int));

  pixbl = realloc(pixbl, naxis * sizeof(double));
  pixbr = realloc(pixbr, naxis * sizeof(double));
  pixul = realloc(pixul, naxis * sizeof(double));
  pixur = realloc(pixur, naxis * sizeof(double));
  wbl = realloc(wbl, naxis * sizeof(double));
  wbr = realloc(wbr, naxis * sizeof(double));
  wul = realloc(wul, naxis * sizeof(double));
  wur = realloc(wur, naxis * sizeof(double));

  for (n = 0; n < naxis; n++) {
    pixbl[n] = 1;
  }
  if ((status = wcsp2s(wcsi, 1, naxis, pixbl, imcrd, &phi, &theta,
                               wbl, st))) {
            wcsperr(wcsi, "");
  } else {
    lgmin = wbl[0];
    lgmax = wbl[0];
    lnmin = wbl[1];
    lnmax = wbl[1];    
  }
  for (n = 0; n < naxis; n++) {
    pixbr[n] = 1;
  }
  pixbr[0] = naxes[0];
  if ((status = wcsp2s(wcsi, 1, naxis, pixbr, imcrd, &phi, &theta,
                               wbr, st))) {
            wcsperr(wcsi, "");
  } else {
    if (lgmin > wbr[0]) {
      lgmin = wbr[0];
    }
    if (lgmax < wbr[0]) {
      lgmax = wbr[0];
    }
    if (lnmin > wbr[1]) {
      lnmin = wbr[1];
    }
    if (lnmax < wbr[1]) {
      lnmax = wbl[1];
    }
  }
  for (n = 0; n < naxis; n++) {
    pixul[n] = 1;
  }
  pixul[1] = naxes[1];
  if ((status = wcsp2s(wcsi, 1, naxis, pixul, imcrd, &phi, &theta,
                               wul, st))) {
            wcsperr(wcsi, "");
  } else {
    if (lgmin > wul[0]) {
      lgmin = wul[0];
    }
    if (lgmax < wul[0]) {
      lgmax = wul[0];
    }
    if (lnmin > wul[1]) {
      lnmin = wul[1];
    }
    if (lnmax < wul[1]) {
      lnmax = wul[1];
    }
  }

  for (n = 0; n < naxis; n++) {
    pixur[n] = 1;
  }
  pixur[0] = naxes[0];
  pixur[1] = naxes[1];
  if ((status = wcsp2s(wcsi, 1, naxis, pixur, imcrd, &phi, &theta,
                               wur, st))) {
            wcsperr(wcsi, "");
  } else {
    if (lgmin > wur[0]) {
      lgmin = wur[0];
    }
    if (lgmax < wur[0]) {
      lgmax = wur[0];
    }
    if (lnmin > wur[1]) {
      lnmin = wur[1];
    }
    if (lnmax < wur[1]) {
      lnmax = wur[1];
    }
  }

  /* Get WCSNAME out of the wcsprm struct. */
  //strcpy(wcsnameout, wcsout->wcsname);
  //printf("%s\n",wcsnameout);

  /* Transform coordinates */
  nelem = wcsout->naxis;

  world  = realloc(world,  nelem * sizeof(double));
  imgcrd = realloc(imgcrd, nelem * sizeof(double));
  pixcrd = realloc(pixcrd, nelem * sizeof(double));
  stat   = realloc(stat,   nelem * sizeof(int));
  imgcrdout = realloc(imgcrdout, naxis * sizeof(double));
  pixcrdout = realloc(pixcrdout, naxis * sizeof(double));
  statout   = realloc(statout,   naxis * sizeof(int));

  worldin  = realloc(worldin,  naxis * sizeof(double));

  if (naxis<3) {
     n3 = 1;
  } else {
     n3 = naxes[2];
  }

  if (n1<1) {
    if (wcsout->cdelt[0] != 0 && wcsout->cdelt[1] != 0) {
      res = wcsout->cdelt[0];
      if (wcsout->crpix[0] != 0. && wcsout->crpix[1] != 0.) {
        n1 = (long)(2*wcsout->crpix[0]) - 1; 
        n2 = (long)(2*wcsout->crpix[1]) - 1; 
      } else {
        n1 = (long)((lgmax - lgmin) / res);
        res = wcsout->cdelt[1];
        n2 = (long)((lnmax - lnmin) / res);
        wcsout->crpix[0] = 0.5 + n1/2.;
        wcsout->crpix[1] = 0.5 + n2/2.;
        wcsout->crval[0] = (lgmax - lgmin) /2.;
        wcsout->crval[1] = (lnmax - lnmin) /2.;
      }
    } else {
      fprintf(stderr, "ERROR: ...");
      return 1;      
    }
  }

  /* Create the output file */
  if (fits_create_file(&fptrout, outfile, &status)) {
      fprintf(stderr, "ERROR: cannot create new file \"%s\".\n", outfile);
      return 1;
  } else {
    naxesout = malloc(naxis*sizeof(long));
    naxesout[0] = n1;
    naxesout[1] = n2;
    if (n3>1) {
      naxesout[2] = n3;
    }
    if (fits_create_img(fptrout, bitpix, naxis, naxesout, &status)) {
        fprintf(stderr, "ERROR: cannot create new image in \"%s\".\n", outfile);
        return 1;
    } else {
    }
  }

  for (z=0; z<n3; z++) {
    for (j=0; j<n2; j++) {
      for (i=0; i<n1; i++) {
          pixcrd[0] = (double)(i+1);
          pixcrd[1] = (double)(j+1);
          
          if ((status = wcsp2s(wcsout, 1, nelem, pixcrd, imgcrd, &phi, &theta,
                               world, stat))) {
            //wcsperr(wcsout, "");

          } else {

            //wcsprintf("\nWorld: ");
            for (n = 0; n < nelem; n++) {
            //    wcsprintf("%s%14.9g", n?", ":"", world[n]);
                worldin[n] = world[n];
            }
            if (n3 - 1) {
              worldin[2] = (double)(z+1);
            }
            
            if ((status = wcss2p(wcsi, 1, naxis, worldin, &phi, &theta, imgcrdout,
                               pixcrdout, statout))) {
              //wcsperr(wcsi, "");
            } else {

              wcsprintf("\nInput: ");
              for (n = 0; n < naxis; n++) {
                  wcsprintf("%s%14.9g", n?", ":"", pixcrdout[n]);
              }
              wcsprintf("\n");
            }
          }
      
      }
    }
  }

  fits_close_file(fptrout, &status);
  fits_close_file(fptr, &status);

  /* Defeat spurious reporting of memory leaks. */
  wcsvfree(&nwcs, &wcs);
  wcsvfree(&nwcsout, &wcsout);
  free(world);
  free(worldin);
  free(imgcrd);
  free(imgcrdout);
  free(pixcrd);
  free(pixcrdout);
  free(stat);
  free(statout);
  free(st);

  return 0;

fitserr:
  fits_report_error(stderr, status);
  fits_close_file(fptr, &status);
  return 1;

}
