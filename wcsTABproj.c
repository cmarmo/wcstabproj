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
"                  output image.\n";


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
  int   i, alts[27], status, nkeyrec, hdutype, relax, ctrl, nreject, nwcs,
        ialt, nfound, bitpix, nwcsout, nkeys, nn;
  long  firstpix, ipix, lastpix, naxis, *naxes = 0x0, nside = 0, repeat,
        offset;

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

  /* Get the output WCS structure */
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
  }

  /* Get WCSNAME out of the wcsprm struct. */
  //strcpy(wcsnameout, wcsout->wcsname);
  //printf("%s\n",wcsnameout);

  /* Transform coordinates */
  

  fits_close_file(fptr, &status);

  /* Defeat spurious reporting of memory leaks. */
  wcsvfree(&nwcs, &wcs);
  wcsvfree(&nwcsout, &wcsout);
  //free(world);
  //free(imgcrd);
  //free(pixcrd);
  //free(stat);

  return 0;

fitserr:
  fits_report_error(stderr, status);
  fits_close_file(fptr, &status);
  return 1;

}
