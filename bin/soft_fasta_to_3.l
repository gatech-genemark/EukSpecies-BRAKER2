%{

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char* seqname=NULL;
static long beg=0L;
static long end=0L;
static int in_lower=0;
static void print_bed()
	{
	printf("%s\t%ld\t%ld\n",seqname,beg,end);
	}
%}

%option noyywrap

%%

<<EOF>>  {if(in_lower) print_bed(); return 0;}
^>.*\n {if(in_lower) print_bed(); free(seqname);seqname=strndup(&yytext[1],yyleng-2);beg=end=0;in_lower=0;}
[A-Z]+ {if(in_lower) {print_bed(); in_lower=0;} end+=yyleng;}
[a-z]+ {if(!in_lower) {beg=end;in_lower=1;} end+=yyleng;}
[ \t\n\r]+  ;
.	{fprintf(stderr,"ERROR\n"); exit(-1);}

%%

int main(int argc,char** argv)
	{
	yylex();
	return 0;
	}
