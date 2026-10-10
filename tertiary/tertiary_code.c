#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdint.h>

#define MAX_CHAR 512

int tcode[210][1500];

int hex2bin(char hex, int* bin)
{
	int val;

	if (hex >= '0' && hex <= '9') {
		val = hex - '0';
	}
	else if (hex >= 'A' && hex <= 'F') {
		val = hex - 'A' + 10;
	}
	else if (hex >= 'a' && hex <= 'f') {
		val = hex - 'a' + 10;
	}
	else { // Invalid hex character
		return 0;
	}

	for (int i = 3; i >= 0; i--) {
		bin[i] = val % 2;
		val /= 2;
	}

	return 1;
}

int readTertiary(char* fname)
{
	FILE* fp;
	char str[MAX_CHAR];
	int nsv;

	if (NULL == (fp = fopen(fname, "rt"))) {
		return 0;
	}

	// Ignore the first line
	fgets(str, MAX_CHAR, fp);

	for (nsv = 0; nsv < 210; nsv++) {

		if (NULL == fgets(str, MAX_CHAR, fp))
			break;

		for (int i = 0; i < 375; i++) {
			hex2bin(str[i + 1], tcode[nsv] + i * 4); // Skip the first charactor
		}
	}

	fclose(fp);

	return nsv;
}

int main()
{
	FILE* fp;

	if (!readTertiary("008_Weil1500hex210prns.txt")) {
		printf("ERROR: Failed to read tertiary code file.\n");
		exit(1);
	}

	if (NULL == (fp = fopen("tcode.txt", "wt"))) {
		printf("ERROR: Failed to open output file.\n");
		exit(1);
	}

	for (int i = 0; i < 210; i++) {
		for (int j = 0; j < 1500; j++) {
			fprintf(fp, "%1d", tcode[i][j]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);

	return 0;
}
