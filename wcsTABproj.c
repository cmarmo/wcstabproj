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
"  -a<alt>      Specify an alternate coordinate representation (ignored if\n"
"               there is only one).  Can also be specified as a 0-relative\n"
"               index in the range 0 to 26, where alternates are sequenced\n"
"               alphabetically following the primary representation.\n"
"  -p<code>     Specify the projection for the output image. Must use the \n"
"               standard FITS WCS codes.\n"
"  -r<deg/pix>  Specify the resolution (in deg/pix) for the output image.\n"
"               It will be translated in CDELT keyword.\n";


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

  int   i, alts[27], status, nkeyrec, hdutype, relax, ctrl, nreject, nwcs,
        ialt, nfound;
  long  firstpix, ipix, lastpix, naxis, *naxes = 0x0, nside = 0, repeat;

  char  *infile;	/* Input file.                                */
  char  *outfile;	/* Output file.                               */
  char  *alt = 0x0, *header, *proj, wcsname[72];

  double res;

  struct wcsprm *wcs, *wcsi = 0x0;
  fitsfile *fptr;

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

    case 'p':
      /* Select a projection. */
      proj = argv[i]+2;
      break;

    case 'r':
      /* Select a resolution. */
      res = atof(argv[i]+2);
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

  fits_close_file(fptr, &status);

  /* Defeat spurious reporting of memory leaks. */
  wcsvfree(&nwcs, &wcs);
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
