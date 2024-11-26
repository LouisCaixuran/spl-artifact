/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison implementation for Yacc-like parsers in C

   Copyright (C) 1984, 1989-1990, 2000-2015, 2018-2021 Free Software Foundation,
   Inc.

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* C LALR(1) parser skeleton written by Richard Stallman, by
   simplifying the original so-called "semantic" parser.  */

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

/* All symbols defined below should begin with yy or YY, to avoid
   infringing on user name space.  This should be done even for local
   variables, as they might otherwise be expanded by user macros.
   There are some unavoidable exceptions within include files to
   define necessary library symbols; they are noted "INFRINGES ON
   USER NAME SPACE" below.  */

/* Identify Bison output, and Bison version.  */
#define YYBISON 30802

/* Bison version string.  */
#define YYBISON_VERSION "3.8.2"

/* Skeleton name.  */
#define YYSKELETON_NAME "yacc.c"

/* Pure parsers.  */
#define YYPURE 0

/* Push parsers.  */
#define YYPUSH 0

/* Pull parsers.  */
#define YYPULL 1




/* First part of user prologue.  */
#line 27 "SDCC.y"

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "SDCCglobl.h"
#include "SDCCattr.h"
#include "SDCCsymt.h"
#include "SDCChasht.h"
#include "SDCCval.h"
#include "SDCCmem.h"
#include "SDCCast.h"
#include "port.h"
#include "newalloc.h"
#include "SDCCerr.h"
#include "SDCCutil.h"
#include "SDCCbtree.h"
#include "SDCCopt.h"
#include "dbuf_string.h"

extern int yyerror (char *);
extern FILE     *yyin;
long NestLevel = 0;     /* current NestLevel       */
int stackPtr  = 1;      /* stack pointer           */
int xstackPtr = 0;      /* xstack pointer          */
int reentrant = 0;
int blockNo   = 0;      /* sequential block number  */
int currBlockno=0;
int inCriticalFunction = 0;
int inCriticalBlock = 0;
int seqPointNo= 1;      /* sequence point number */
int ignoreTypedefType=0;
extern int yylex();
int yyparse(void);
extern int noLineno;
char lbuff[1024];       /* local buffer */
char function_name[256] = {0};

/* break & continue stacks */
STACK_DCL(continueStack  ,symbol *,MAX_NEST_LEVEL)
STACK_DCL(breakStack  ,symbol *,MAX_NEST_LEVEL)
STACK_DCL(forStack  ,symbol *,MAX_NEST_LEVEL)
STACK_DCL(swStk   ,ast   *,MAX_NEST_LEVEL)
STACK_DCL(blockNum,int,MAX_NEST_LEVEL*3)

value *cenum = NULL;        /* current enumeration  type chain*/
bool uselessDecl = true;

#define YYDEBUG 1


#line 122 "y.tab.c"

# ifndef YY_CAST
#  ifdef __cplusplus
#   define YY_CAST(Type, Val) static_cast<Type> (Val)
#   define YY_REINTERPRET_CAST(Type, Val) reinterpret_cast<Type> (Val)
#  else
#   define YY_CAST(Type, Val) ((Type) (Val))
#   define YY_REINTERPRET_CAST(Type, Val) ((Type) (Val))
#  endif
# endif
# ifndef YY_NULLPTR
#  if defined __cplusplus
#   if 201103L <= __cplusplus
#    define YY_NULLPTR nullptr
#   else
#    define YY_NULLPTR 0
#   endif
#  else
#   define YY_NULLPTR ((void*)0)
#  endif
# endif

/* Use api.header.include to #include this header
   instead of duplicating it here.  */
#ifndef YY_YY_Y_TAB_H_INCLUDED
# define YY_YY_Y_TAB_H_INCLUDED
/* Debug traces.  */
#ifndef YYDEBUG
# define YYDEBUG 0
#endif
#if YYDEBUG
extern int yydebug;
#endif

/* Token kinds.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
  enum yytokentype
  {
    YYEMPTY = -2,
    YYEOF = 0,                     /* "end of file"  */
    YYerror = 256,                 /* error  */
    YYUNDEF = 257,                 /* "invalid token"  */
    IDENTIFIER = 258,              /* IDENTIFIER  */
    TYPE_NAME = 259,               /* TYPE_NAME  */
    ADDRSPACE_NAME = 260,          /* ADDRSPACE_NAME  */
    CONSTANT = 261,                /* CONSTANT  */
    SIZEOF = 262,                  /* SIZEOF  */
    OFFSETOF = 263,                /* OFFSETOF  */
    PTR_OP = 264,                  /* PTR_OP  */
    INC_OP = 265,                  /* INC_OP  */
    DEC_OP = 266,                  /* DEC_OP  */
    LEFT_OP = 267,                 /* LEFT_OP  */
    RIGHT_OP = 268,                /* RIGHT_OP  */
    LE_OP = 269,                   /* LE_OP  */
    GE_OP = 270,                   /* GE_OP  */
    EQ_OP = 271,                   /* EQ_OP  */
    NE_OP = 272,                   /* NE_OP  */
    AND_OP = 273,                  /* AND_OP  */
    OR_OP = 274,                   /* OR_OP  */
    ATTR_START = 275,              /* ATTR_START  */
    TOK_SEP = 276,                 /* TOK_SEP  */
    MUL_ASSIGN = 277,              /* MUL_ASSIGN  */
    DIV_ASSIGN = 278,              /* DIV_ASSIGN  */
    MOD_ASSIGN = 279,              /* MOD_ASSIGN  */
    ADD_ASSIGN = 280,              /* ADD_ASSIGN  */
    SUB_ASSIGN = 281,              /* SUB_ASSIGN  */
    LEFT_ASSIGN = 282,             /* LEFT_ASSIGN  */
    RIGHT_ASSIGN = 283,            /* RIGHT_ASSIGN  */
    AND_ASSIGN = 284,              /* AND_ASSIGN  */
    XOR_ASSIGN = 285,              /* XOR_ASSIGN  */
    OR_ASSIGN = 286,               /* OR_ASSIGN  */
    TYPEDEF = 287,                 /* TYPEDEF  */
    EXTERN = 288,                  /* EXTERN  */
    STATIC = 289,                  /* STATIC  */
    AUTO = 290,                    /* AUTO  */
    REGISTER = 291,                /* REGISTER  */
    CONSTEXPR = 292,               /* CONSTEXPR  */
    CODE = 293,                    /* CODE  */
    EEPROM = 294,                  /* EEPROM  */
    INTERRUPT = 295,               /* INTERRUPT  */
    SFR = 296,                     /* SFR  */
    SFR16 = 297,                   /* SFR16  */
    SFR32 = 298,                   /* SFR32  */
    ADDRESSMOD = 299,              /* ADDRESSMOD  */
    AT = 300,                      /* AT  */
    SBIT = 301,                    /* SBIT  */
    REENTRANT = 302,               /* REENTRANT  */
    USING = 303,                   /* USING  */
    XDATA = 304,                   /* XDATA  */
    DATA = 305,                    /* DATA  */
    IDATA = 306,                   /* IDATA  */
    PDATA = 307,                   /* PDATA  */
    ELLIPSIS = 308,                /* ELLIPSIS  */
    CRITICAL = 309,                /* CRITICAL  */
    NONBANKED = 310,               /* NONBANKED  */
    BANKED = 311,                  /* BANKED  */
    SHADOWREGS = 312,              /* SHADOWREGS  */
    SD_WPARAM = 313,               /* SD_WPARAM  */
    SD_BOOL = 314,                 /* SD_BOOL  */
    SD_CHAR = 315,                 /* SD_CHAR  */
    SD_SHORT = 316,                /* SD_SHORT  */
    SD_INT = 317,                  /* SD_INT  */
    SD_LONG = 318,                 /* SD_LONG  */
    SIGNED = 319,                  /* SIGNED  */
    UNSIGNED = 320,                /* UNSIGNED  */
    SD_FLOAT = 321,                /* SD_FLOAT  */
    DOUBLE = 322,                  /* DOUBLE  */
    FIXED16X16 = 323,              /* FIXED16X16  */
    SD_CONST = 324,                /* SD_CONST  */
    VOLATILE = 325,                /* VOLATILE  */
    SD_VOID = 326,                 /* SD_VOID  */
    BIT = 327,                     /* BIT  */
    COMPLEX = 328,                 /* COMPLEX  */
    IMAGINARY = 329,               /* IMAGINARY  */
    STRUCT = 330,                  /* STRUCT  */
    UNION = 331,                   /* UNION  */
    ENUM = 332,                    /* ENUM  */
    RANGE = 333,                   /* RANGE  */
    SD_FAR = 334,                  /* SD_FAR  */
    CASE = 335,                    /* CASE  */
    DEFAULT = 336,                 /* DEFAULT  */
    IF = 337,                      /* IF  */
    ELSE = 338,                    /* ELSE  */
    SWITCH = 339,                  /* SWITCH  */
    WHILE = 340,                   /* WHILE  */
    DO = 341,                      /* DO  */
    FOR = 342,                     /* FOR  */
    GOTO = 343,                    /* GOTO  */
    CONTINUE = 344,                /* CONTINUE  */
    BREAK = 345,                   /* BREAK  */
    RETURN = 346,                  /* RETURN  */
    NAKED = 347,                   /* NAKED  */
    JAVANATIVE = 348,              /* JAVANATIVE  */
    OVERLAY = 349,                 /* OVERLAY  */
    TRAP = 350,                    /* TRAP  */
    STRING_LITERAL = 351,          /* STRING_LITERAL  */
    INLINEASM = 352,               /* INLINEASM  */
    FUNC = 353,                    /* FUNC  */
    IFX = 354,                     /* IFX  */
    ADDRESS_OF = 355,              /* ADDRESS_OF  */
    GET_VALUE_AT_ADDRESS = 356,    /* GET_VALUE_AT_ADDRESS  */
    SET_VALUE_AT_ADDRESS = 357,    /* SET_VALUE_AT_ADDRESS  */
    SPIL = 358,                    /* SPIL  */
    UNSPIL = 359,                  /* UNSPIL  */
    GETABIT = 360,                 /* GETABIT  */
    GETBYTE = 361,                 /* GETBYTE  */
    GETWORD = 362,                 /* GETWORD  */
    BITWISEAND = 363,              /* BITWISEAND  */
    UNARYMINUS = 364,              /* UNARYMINUS  */
    IPUSH = 365,                   /* IPUSH  */
    IPUSH_VALUE_AT_ADDRESS = 366,  /* IPUSH_VALUE_AT_ADDRESS  */
    IPOP = 367,                    /* IPOP  */
    PCALL = 368,                   /* PCALL  */
    ENDFUNCTION = 369,             /* ENDFUNCTION  */
    JUMPTABLE = 370,               /* JUMPTABLE  */
    RRC = 371,                     /* RRC  */
    RLC = 372,                     /* RLC  */
    CAST = 373,                    /* CAST  */
    CALL = 374,                    /* CALL  */
    PARAM = 375,                   /* PARAM  */
    NULLOP = 376,                  /* NULLOP  */
    BLOCK = 377,                   /* BLOCK  */
    LABEL = 378,                   /* LABEL  */
    RECEIVE = 379,                 /* RECEIVE  */
    SEND = 380,                    /* SEND  */
    ARRAYINIT = 381,               /* ARRAYINIT  */
    DUMMY_READ_VOLATILE = 382,     /* DUMMY_READ_VOLATILE  */
    ENDCRITICAL = 383,             /* ENDCRITICAL  */
    SWAP = 384,                    /* SWAP  */
    INLINE = 385,                  /* INLINE  */
    RESTRICT = 386,                /* RESTRICT  */
    SMALLC = 387,                  /* SMALLC  */
    RAISONANCE = 388,              /* RAISONANCE  */
    IAR = 389,                     /* IAR  */
    COSMIC = 390,                  /* COSMIC  */
    SDCCCALL = 391,                /* SDCCCALL  */
    PRESERVES_REGS = 392,          /* PRESERVES_REGS  */
    Z88DK_FASTCALL = 393,          /* Z88DK_FASTCALL  */
    Z88DK_CALLEE = 394,            /* Z88DK_CALLEE  */
    Z88DK_SHORTCALL = 395,         /* Z88DK_SHORTCALL  */
    Z88DK_PARAMS_OFFSET = 396,     /* Z88DK_PARAMS_OFFSET  */
    ALIGNAS = 397,                 /* ALIGNAS  */
    ALIGNOF = 398,                 /* ALIGNOF  */
    ATOMIC = 399,                  /* ATOMIC  */
    GENERIC = 400,                 /* GENERIC  */
    NORETURN = 401,                /* NORETURN  */
    STATIC_ASSERT = 402,           /* STATIC_ASSERT  */
    THREAD_LOCAL = 403,            /* THREAD_LOCAL  */
    TOKEN_FALSE = 404,             /* TOKEN_FALSE  */
    TOKEN_TRUE = 405,              /* TOKEN_TRUE  */
    NULLPTR = 406,                 /* NULLPTR  */
    TYPEOF = 407,                  /* TYPEOF  */
    TYPEOF_UNQUAL = 408,           /* TYPEOF_UNQUAL  */
    SD_BITINT = 409,               /* SD_BITINT  */
    DECIMAL32 = 410,               /* DECIMAL32  */
    DECIMAL64 = 411,               /* DECIMAL64  */
    DECIMAL128 = 412,              /* DECIMAL128  */
    ASM = 413,                     /* ASM  */
    GENERIC_ASSOC_LIST = 414,      /* GENERIC_ASSOC_LIST  */
    GENERIC_ASSOCIATION = 415      /* GENERIC_ASSOCIATION  */
  };
  typedef enum yytokentype yytoken_kind_t;
#endif
/* Token kinds.  */
#define YYEMPTY -2
#define YYEOF 0
#define YYerror 256
#define YYUNDEF 257
#define IDENTIFIER 258
#define TYPE_NAME 259
#define ADDRSPACE_NAME 260
#define CONSTANT 261
#define SIZEOF 262
#define OFFSETOF 263
#define PTR_OP 264
#define INC_OP 265
#define DEC_OP 266
#define LEFT_OP 267
#define RIGHT_OP 268
#define LE_OP 269
#define GE_OP 270
#define EQ_OP 271
#define NE_OP 272
#define AND_OP 273
#define OR_OP 274
#define ATTR_START 275
#define TOK_SEP 276
#define MUL_ASSIGN 277
#define DIV_ASSIGN 278
#define MOD_ASSIGN 279
#define ADD_ASSIGN 280
#define SUB_ASSIGN 281
#define LEFT_ASSIGN 282
#define RIGHT_ASSIGN 283
#define AND_ASSIGN 284
#define XOR_ASSIGN 285
#define OR_ASSIGN 286
#define TYPEDEF 287
#define EXTERN 288
#define STATIC 289
#define AUTO 290
#define REGISTER 291
#define CONSTEXPR 292
#define CODE 293
#define EEPROM 294
#define INTERRUPT 295
#define SFR 296
#define SFR16 297
#define SFR32 298
#define ADDRESSMOD 299
#define AT 300
#define SBIT 301
#define REENTRANT 302
#define USING 303
#define XDATA 304
#define DATA 305
#define IDATA 306
#define PDATA 307
#define ELLIPSIS 308
#define CRITICAL 309
#define NONBANKED 310
#define BANKED 311
#define SHADOWREGS 312
#define SD_WPARAM 313
#define SD_BOOL 314
#define SD_CHAR 315
#define SD_SHORT 316
#define SD_INT 317
#define SD_LONG 318
#define SIGNED 319
#define UNSIGNED 320
#define SD_FLOAT 321
#define DOUBLE 322
#define FIXED16X16 323
#define SD_CONST 324
#define VOLATILE 325
#define SD_VOID 326
#define BIT 327
#define COMPLEX 328
#define IMAGINARY 329
#define STRUCT 330
#define UNION 331
#define ENUM 332
#define RANGE 333
#define SD_FAR 334
#define CASE 335
#define DEFAULT 336
#define IF 337
#define ELSE 338
#define SWITCH 339
#define WHILE 340
#define DO 341
#define FOR 342
#define GOTO 343
#define CONTINUE 344
#define BREAK 345
#define RETURN 346
#define NAKED 347
#define JAVANATIVE 348
#define OVERLAY 349
#define TRAP 350
#define STRING_LITERAL 351
#define INLINEASM 352
#define FUNC 353
#define IFX 354
#define ADDRESS_OF 355
#define GET_VALUE_AT_ADDRESS 356
#define SET_VALUE_AT_ADDRESS 357
#define SPIL 358
#define UNSPIL 359
#define GETABIT 360
#define GETBYTE 361
#define GETWORD 362
#define BITWISEAND 363
#define UNARYMINUS 364
#define IPUSH 365
#define IPUSH_VALUE_AT_ADDRESS 366
#define IPOP 367
#define PCALL 368
#define ENDFUNCTION 369
#define JUMPTABLE 370
#define RRC 371
#define RLC 372
#define CAST 373
#define CALL 374
#define PARAM 375
#define NULLOP 376
#define BLOCK 377
#define LABEL 378
#define RECEIVE 379
#define SEND 380
#define ARRAYINIT 381
#define DUMMY_READ_VOLATILE 382
#define ENDCRITICAL 383
#define SWAP 384
#define INLINE 385
#define RESTRICT 386
#define SMALLC 387
#define RAISONANCE 388
#define IAR 389
#define COSMIC 390
#define SDCCCALL 391
#define PRESERVES_REGS 392
#define Z88DK_FASTCALL 393
#define Z88DK_CALLEE 394
#define Z88DK_SHORTCALL 395
#define Z88DK_PARAMS_OFFSET 396
#define ALIGNAS 397
#define ALIGNOF 398
#define ATOMIC 399
#define GENERIC 400
#define NORETURN 401
#define STATIC_ASSERT 402
#define THREAD_LOCAL 403
#define TOKEN_FALSE 404
#define TOKEN_TRUE 405
#define NULLPTR 406
#define TYPEOF 407
#define TYPEOF_UNQUAL 408
#define SD_BITINT 409
#define DECIMAL32 410
#define DECIMAL64 411
#define DECIMAL128 412
#define ASM 413
#define GENERIC_ASSOC_LIST 414
#define GENERIC_ASSOCIATION 415

/* Value type.  */
#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
union YYSTYPE
{
#line 79 "SDCC.y"

    attribute  *attr;       /* attribute                              */
    symbol     *sym;        /* symbol table pointer                   */
    structdef  *sdef;       /* structure definition                   */
    char       yychar[SDCC_NAME_MAX+1];
    sym_link   *lnk;        /* declarator  or specifier               */
    int        yyint;       /* integer value returned                 */
    value      *val;        /* for integer constant                   */
    initList   *ilist;      /* initial list                           */
    designation*dsgn;       /* designator                             */
    const char *yystr;      /* pointer to dynamicaly allocated string */
    ast        *asts;       /* expression tree                        */

#line 509 "y.tab.c"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;


int yyparse (void);


#endif /* !YY_YY_Y_TAB_H_INCLUDED  */
/* Symbol kind.  */
enum yysymbol_kind_t
{
  YYSYMBOL_YYEMPTY = -2,
  YYSYMBOL_YYEOF = 0,                      /* "end of file"  */
  YYSYMBOL_YYerror = 1,                    /* error  */
  YYSYMBOL_YYUNDEF = 2,                    /* "invalid token"  */
  YYSYMBOL_IDENTIFIER = 3,                 /* IDENTIFIER  */
  YYSYMBOL_TYPE_NAME = 4,                  /* TYPE_NAME  */
  YYSYMBOL_ADDRSPACE_NAME = 5,             /* ADDRSPACE_NAME  */
  YYSYMBOL_CONSTANT = 6,                   /* CONSTANT  */
  YYSYMBOL_SIZEOF = 7,                     /* SIZEOF  */
  YYSYMBOL_OFFSETOF = 8,                   /* OFFSETOF  */
  YYSYMBOL_PTR_OP = 9,                     /* PTR_OP  */
  YYSYMBOL_INC_OP = 10,                    /* INC_OP  */
  YYSYMBOL_DEC_OP = 11,                    /* DEC_OP  */
  YYSYMBOL_LEFT_OP = 12,                   /* LEFT_OP  */
  YYSYMBOL_RIGHT_OP = 13,                  /* RIGHT_OP  */
  YYSYMBOL_LE_OP = 14,                     /* LE_OP  */
  YYSYMBOL_GE_OP = 15,                     /* GE_OP  */
  YYSYMBOL_EQ_OP = 16,                     /* EQ_OP  */
  YYSYMBOL_NE_OP = 17,                     /* NE_OP  */
  YYSYMBOL_AND_OP = 18,                    /* AND_OP  */
  YYSYMBOL_OR_OP = 19,                     /* OR_OP  */
  YYSYMBOL_ATTR_START = 20,                /* ATTR_START  */
  YYSYMBOL_TOK_SEP = 21,                   /* TOK_SEP  */
  YYSYMBOL_MUL_ASSIGN = 22,                /* MUL_ASSIGN  */
  YYSYMBOL_DIV_ASSIGN = 23,                /* DIV_ASSIGN  */
  YYSYMBOL_MOD_ASSIGN = 24,                /* MOD_ASSIGN  */
  YYSYMBOL_ADD_ASSIGN = 25,                /* ADD_ASSIGN  */
  YYSYMBOL_SUB_ASSIGN = 26,                /* SUB_ASSIGN  */
  YYSYMBOL_LEFT_ASSIGN = 27,               /* LEFT_ASSIGN  */
  YYSYMBOL_RIGHT_ASSIGN = 28,              /* RIGHT_ASSIGN  */
  YYSYMBOL_AND_ASSIGN = 29,                /* AND_ASSIGN  */
  YYSYMBOL_XOR_ASSIGN = 30,                /* XOR_ASSIGN  */
  YYSYMBOL_OR_ASSIGN = 31,                 /* OR_ASSIGN  */
  YYSYMBOL_TYPEDEF = 32,                   /* TYPEDEF  */
  YYSYMBOL_EXTERN = 33,                    /* EXTERN  */
  YYSYMBOL_STATIC = 34,                    /* STATIC  */
  YYSYMBOL_AUTO = 35,                      /* AUTO  */
  YYSYMBOL_REGISTER = 36,                  /* REGISTER  */
  YYSYMBOL_CONSTEXPR = 37,                 /* CONSTEXPR  */
  YYSYMBOL_CODE = 38,                      /* CODE  */
  YYSYMBOL_EEPROM = 39,                    /* EEPROM  */
  YYSYMBOL_INTERRUPT = 40,                 /* INTERRUPT  */
  YYSYMBOL_SFR = 41,                       /* SFR  */
  YYSYMBOL_SFR16 = 42,                     /* SFR16  */
  YYSYMBOL_SFR32 = 43,                     /* SFR32  */
  YYSYMBOL_ADDRESSMOD = 44,                /* ADDRESSMOD  */
  YYSYMBOL_AT = 45,                        /* AT  */
  YYSYMBOL_SBIT = 46,                      /* SBIT  */
  YYSYMBOL_REENTRANT = 47,                 /* REENTRANT  */
  YYSYMBOL_USING = 48,                     /* USING  */
  YYSYMBOL_XDATA = 49,                     /* XDATA  */
  YYSYMBOL_DATA = 50,                      /* DATA  */
  YYSYMBOL_IDATA = 51,                     /* IDATA  */
  YYSYMBOL_PDATA = 52,                     /* PDATA  */
  YYSYMBOL_ELLIPSIS = 53,                  /* ELLIPSIS  */
  YYSYMBOL_CRITICAL = 54,                  /* CRITICAL  */
  YYSYMBOL_NONBANKED = 55,                 /* NONBANKED  */
  YYSYMBOL_BANKED = 56,                    /* BANKED  */
  YYSYMBOL_SHADOWREGS = 57,                /* SHADOWREGS  */
  YYSYMBOL_SD_WPARAM = 58,                 /* SD_WPARAM  */
  YYSYMBOL_SD_BOOL = 59,                   /* SD_BOOL  */
  YYSYMBOL_SD_CHAR = 60,                   /* SD_CHAR  */
  YYSYMBOL_SD_SHORT = 61,                  /* SD_SHORT  */
  YYSYMBOL_SD_INT = 62,                    /* SD_INT  */
  YYSYMBOL_SD_LONG = 63,                   /* SD_LONG  */
  YYSYMBOL_SIGNED = 64,                    /* SIGNED  */
  YYSYMBOL_UNSIGNED = 65,                  /* UNSIGNED  */
  YYSYMBOL_SD_FLOAT = 66,                  /* SD_FLOAT  */
  YYSYMBOL_DOUBLE = 67,                    /* DOUBLE  */
  YYSYMBOL_FIXED16X16 = 68,                /* FIXED16X16  */
  YYSYMBOL_SD_CONST = 69,                  /* SD_CONST  */
  YYSYMBOL_VOLATILE = 70,                  /* VOLATILE  */
  YYSYMBOL_SD_VOID = 71,                   /* SD_VOID  */
  YYSYMBOL_BIT = 72,                       /* BIT  */
  YYSYMBOL_COMPLEX = 73,                   /* COMPLEX  */
  YYSYMBOL_IMAGINARY = 74,                 /* IMAGINARY  */
  YYSYMBOL_STRUCT = 75,                    /* STRUCT  */
  YYSYMBOL_UNION = 76,                     /* UNION  */
  YYSYMBOL_ENUM = 77,                      /* ENUM  */
  YYSYMBOL_RANGE = 78,                     /* RANGE  */
  YYSYMBOL_SD_FAR = 79,                    /* SD_FAR  */
  YYSYMBOL_CASE = 80,                      /* CASE  */
  YYSYMBOL_DEFAULT = 81,                   /* DEFAULT  */
  YYSYMBOL_IF = 82,                        /* IF  */
  YYSYMBOL_ELSE = 83,                      /* ELSE  */
  YYSYMBOL_SWITCH = 84,                    /* SWITCH  */
  YYSYMBOL_WHILE = 85,                     /* WHILE  */
  YYSYMBOL_DO = 86,                        /* DO  */
  YYSYMBOL_FOR = 87,                       /* FOR  */
  YYSYMBOL_GOTO = 88,                      /* GOTO  */
  YYSYMBOL_CONTINUE = 89,                  /* CONTINUE  */
  YYSYMBOL_BREAK = 90,                     /* BREAK  */
  YYSYMBOL_RETURN = 91,                    /* RETURN  */
  YYSYMBOL_NAKED = 92,                     /* NAKED  */
  YYSYMBOL_JAVANATIVE = 93,                /* JAVANATIVE  */
  YYSYMBOL_OVERLAY = 94,                   /* OVERLAY  */
  YYSYMBOL_TRAP = 95,                      /* TRAP  */
  YYSYMBOL_STRING_LITERAL = 96,            /* STRING_LITERAL  */
  YYSYMBOL_INLINEASM = 97,                 /* INLINEASM  */
  YYSYMBOL_FUNC = 98,                      /* FUNC  */
  YYSYMBOL_IFX = 99,                       /* IFX  */
  YYSYMBOL_ADDRESS_OF = 100,               /* ADDRESS_OF  */
  YYSYMBOL_GET_VALUE_AT_ADDRESS = 101,     /* GET_VALUE_AT_ADDRESS  */
  YYSYMBOL_SET_VALUE_AT_ADDRESS = 102,     /* SET_VALUE_AT_ADDRESS  */
  YYSYMBOL_SPIL = 103,                     /* SPIL  */
  YYSYMBOL_UNSPIL = 104,                   /* UNSPIL  */
  YYSYMBOL_GETABIT = 105,                  /* GETABIT  */
  YYSYMBOL_GETBYTE = 106,                  /* GETBYTE  */
  YYSYMBOL_GETWORD = 107,                  /* GETWORD  */
  YYSYMBOL_BITWISEAND = 108,               /* BITWISEAND  */
  YYSYMBOL_UNARYMINUS = 109,               /* UNARYMINUS  */
  YYSYMBOL_IPUSH = 110,                    /* IPUSH  */
  YYSYMBOL_IPUSH_VALUE_AT_ADDRESS = 111,   /* IPUSH_VALUE_AT_ADDRESS  */
  YYSYMBOL_IPOP = 112,                     /* IPOP  */
  YYSYMBOL_PCALL = 113,                    /* PCALL  */
  YYSYMBOL_ENDFUNCTION = 114,              /* ENDFUNCTION  */
  YYSYMBOL_JUMPTABLE = 115,                /* JUMPTABLE  */
  YYSYMBOL_RRC = 116,                      /* RRC  */
  YYSYMBOL_RLC = 117,                      /* RLC  */
  YYSYMBOL_CAST = 118,                     /* CAST  */
  YYSYMBOL_CALL = 119,                     /* CALL  */
  YYSYMBOL_PARAM = 120,                    /* PARAM  */
  YYSYMBOL_NULLOP = 121,                   /* NULLOP  */
  YYSYMBOL_BLOCK = 122,                    /* BLOCK  */
  YYSYMBOL_LABEL = 123,                    /* LABEL  */
  YYSYMBOL_RECEIVE = 124,                  /* RECEIVE  */
  YYSYMBOL_SEND = 125,                     /* SEND  */
  YYSYMBOL_ARRAYINIT = 126,                /* ARRAYINIT  */
  YYSYMBOL_DUMMY_READ_VOLATILE = 127,      /* DUMMY_READ_VOLATILE  */
  YYSYMBOL_ENDCRITICAL = 128,              /* ENDCRITICAL  */
  YYSYMBOL_SWAP = 129,                     /* SWAP  */
  YYSYMBOL_INLINE = 130,                   /* INLINE  */
  YYSYMBOL_RESTRICT = 131,                 /* RESTRICT  */
  YYSYMBOL_SMALLC = 132,                   /* SMALLC  */
  YYSYMBOL_RAISONANCE = 133,               /* RAISONANCE  */
  YYSYMBOL_IAR = 134,                      /* IAR  */
  YYSYMBOL_COSMIC = 135,                   /* COSMIC  */
  YYSYMBOL_SDCCCALL = 136,                 /* SDCCCALL  */
  YYSYMBOL_PRESERVES_REGS = 137,           /* PRESERVES_REGS  */
  YYSYMBOL_Z88DK_FASTCALL = 138,           /* Z88DK_FASTCALL  */
  YYSYMBOL_Z88DK_CALLEE = 139,             /* Z88DK_CALLEE  */
  YYSYMBOL_Z88DK_SHORTCALL = 140,          /* Z88DK_SHORTCALL  */
  YYSYMBOL_Z88DK_PARAMS_OFFSET = 141,      /* Z88DK_PARAMS_OFFSET  */
  YYSYMBOL_ALIGNAS = 142,                  /* ALIGNAS  */
  YYSYMBOL_ALIGNOF = 143,                  /* ALIGNOF  */
  YYSYMBOL_ATOMIC = 144,                   /* ATOMIC  */
  YYSYMBOL_GENERIC = 145,                  /* GENERIC  */
  YYSYMBOL_NORETURN = 146,                 /* NORETURN  */
  YYSYMBOL_STATIC_ASSERT = 147,            /* STATIC_ASSERT  */
  YYSYMBOL_THREAD_LOCAL = 148,             /* THREAD_LOCAL  */
  YYSYMBOL_TOKEN_FALSE = 149,              /* TOKEN_FALSE  */
  YYSYMBOL_TOKEN_TRUE = 150,               /* TOKEN_TRUE  */
  YYSYMBOL_NULLPTR = 151,                  /* NULLPTR  */
  YYSYMBOL_TYPEOF = 152,                   /* TYPEOF  */
  YYSYMBOL_TYPEOF_UNQUAL = 153,            /* TYPEOF_UNQUAL  */
  YYSYMBOL_SD_BITINT = 154,                /* SD_BITINT  */
  YYSYMBOL_DECIMAL32 = 155,                /* DECIMAL32  */
  YYSYMBOL_DECIMAL64 = 156,                /* DECIMAL64  */
  YYSYMBOL_DECIMAL128 = 157,               /* DECIMAL128  */
  YYSYMBOL_ASM = 158,                      /* ASM  */
  YYSYMBOL_GENERIC_ASSOC_LIST = 159,       /* GENERIC_ASSOC_LIST  */
  YYSYMBOL_GENERIC_ASSOCIATION = 160,      /* GENERIC_ASSOCIATION  */
  YYSYMBOL_161_ = 161,                     /* '('  */
  YYSYMBOL_162_ = 162,                     /* ')'  */
  YYSYMBOL_163_ = 163,                     /* ','  */
  YYSYMBOL_164_ = 164,                     /* ':'  */
  YYSYMBOL_165_ = 165,                     /* '['  */
  YYSYMBOL_166_ = 166,                     /* ']'  */
  YYSYMBOL_167_ = 167,                     /* '.'  */
  YYSYMBOL_168_ = 168,                     /* '{'  */
  YYSYMBOL_169_ = 169,                     /* '}'  */
  YYSYMBOL_170_ = 170,                     /* '&'  */
  YYSYMBOL_171_ = 171,                     /* '*'  */
  YYSYMBOL_172_ = 172,                     /* '+'  */
  YYSYMBOL_173_ = 173,                     /* '-'  */
  YYSYMBOL_174_ = 174,                     /* '~'  */
  YYSYMBOL_175_ = 175,                     /* '!'  */
  YYSYMBOL_176_ = 176,                     /* '/'  */
  YYSYMBOL_177_ = 177,                     /* '%'  */
  YYSYMBOL_178_ = 178,                     /* '<'  */
  YYSYMBOL_179_ = 179,                     /* '>'  */
  YYSYMBOL_180_ = 180,                     /* '^'  */
  YYSYMBOL_181_ = 181,                     /* '|'  */
  YYSYMBOL_182_ = 182,                     /* '?'  */
  YYSYMBOL_183_ = 183,                     /* '='  */
  YYSYMBOL_184_ = 184,                     /* ';'  */
  YYSYMBOL_YYACCEPT = 185,                 /* $accept  */
  YYSYMBOL_primary_expression = 186,       /* primary_expression  */
  YYSYMBOL_predefined_constant = 187,      /* predefined_constant  */
  YYSYMBOL_generic_selection = 188,        /* generic_selection  */
  YYSYMBOL_generic_assoc_list = 189,       /* generic_assoc_list  */
  YYSYMBOL_generic_association = 190,      /* generic_association  */
  YYSYMBOL_postfix_expression = 191,       /* postfix_expression  */
  YYSYMBOL_192_1 = 192,                    /* $@1  */
  YYSYMBOL_193_2 = 193,                    /* $@2  */
  YYSYMBOL_argument_expr_list = 194,       /* argument_expr_list  */
  YYSYMBOL_unary_expression = 195,         /* unary_expression  */
  YYSYMBOL_unary_operator = 196,           /* unary_operator  */
  YYSYMBOL_cast_expression = 197,          /* cast_expression  */
  YYSYMBOL_multiplicative_expression = 198, /* multiplicative_expression  */
  YYSYMBOL_additive_expression = 199,      /* additive_expression  */
  YYSYMBOL_shift_expression = 200,         /* shift_expression  */
  YYSYMBOL_relational_expression = 201,    /* relational_expression  */
  YYSYMBOL_equality_expression = 202,      /* equality_expression  */
  YYSYMBOL_and_expression = 203,           /* and_expression  */
  YYSYMBOL_exclusive_or_expression = 204,  /* exclusive_or_expression  */
  YYSYMBOL_inclusive_or_expression = 205,  /* inclusive_or_expression  */
  YYSYMBOL_logical_and_expr = 206,         /* logical_and_expr  */
  YYSYMBOL_207_3 = 207,                    /* $@3  */
  YYSYMBOL_logical_or_expr = 208,          /* logical_or_expr  */
  YYSYMBOL_209_4 = 209,                    /* $@4  */
  YYSYMBOL_conditional_expr = 210,         /* conditional_expr  */
  YYSYMBOL_211_5 = 211,                    /* $@5  */
  YYSYMBOL_assignment_expr = 212,          /* assignment_expr  */
  YYSYMBOL_assignment_operator = 213,      /* assignment_operator  */
  YYSYMBOL_expression = 214,               /* expression  */
  YYSYMBOL_215_6 = 215,                    /* $@6  */
  YYSYMBOL_expression_opt = 216,           /* expression_opt  */
  YYSYMBOL_constant_expr = 217,            /* constant_expr  */
  YYSYMBOL_declaration = 218,              /* declaration  */
  YYSYMBOL_declaration_specifiers = 219,   /* declaration_specifiers  */
  YYSYMBOL_declaration_specifiers_ = 220,  /* declaration_specifiers_  */
  YYSYMBOL_init_declarator_list = 221,     /* init_declarator_list  */
  YYSYMBOL_init_declarator = 222,          /* init_declarator  */
  YYSYMBOL_attribute_declaration = 223,    /* attribute_declaration  */
  YYSYMBOL_storage_class_specifier = 224,  /* storage_class_specifier  */
  YYSYMBOL_type_specifier = 225,           /* type_specifier  */
  YYSYMBOL_typeof_specifier = 226,         /* typeof_specifier  */
  YYSYMBOL_struct_or_union_specifier = 227, /* struct_or_union_specifier  */
  YYSYMBOL_228_7 = 228,                    /* $@7  */
  YYSYMBOL_struct_or_union = 229,          /* struct_or_union  */
  YYSYMBOL_member_declaration_list = 230,  /* member_declaration_list  */
  YYSYMBOL_member_declaration = 231,       /* member_declaration  */
  YYSYMBOL_type_specifier_qualifier = 232, /* type_specifier_qualifier  */
  YYSYMBOL_member_declarator_list = 233,   /* member_declarator_list  */
  YYSYMBOL_member_declarator = 234,        /* member_declarator  */
  YYSYMBOL_enum_specifier = 235,           /* enum_specifier  */
  YYSYMBOL_enum_comma_opt = 236,           /* enum_comma_opt  */
  YYSYMBOL_enumerator_list = 237,          /* enumerator_list  */
  YYSYMBOL_enumerator = 238,               /* enumerator  */
  YYSYMBOL_type_qualifier = 239,           /* type_qualifier  */
  YYSYMBOL_function_specifier = 240,       /* function_specifier  */
  YYSYMBOL_alignment_specifier = 241,      /* alignment_specifier  */
  YYSYMBOL_declarator = 242,               /* declarator  */
  YYSYMBOL_direct_declarator = 243,        /* direct_declarator  */
  YYSYMBOL_declarator2 = 244,              /* declarator2  */
  YYSYMBOL_array_declarator = 245,         /* array_declarator  */
  YYSYMBOL_declarator2_function_attributes = 246, /* declarator2_function_attributes  */
  YYSYMBOL_function_declarator = 247,      /* function_declarator  */
  YYSYMBOL_248_8 = 248,                    /* $@8  */
  YYSYMBOL_pointer = 249,                  /* pointer  */
  YYSYMBOL_unqualified_pointer = 250,      /* unqualified_pointer  */
  YYSYMBOL_type_qualifier_list = 251,      /* type_qualifier_list  */
  YYSYMBOL_type_qualifier_list_opt = 252,  /* type_qualifier_list_opt  */
  YYSYMBOL_parameter_type_list = 253,      /* parameter_type_list  */
  YYSYMBOL_parameter_list = 254,           /* parameter_list  */
  YYSYMBOL_parameter_declaration = 255,    /* parameter_declaration  */
  YYSYMBOL_abstract_declarator = 256,      /* abstract_declarator  */
  YYSYMBOL_direct_abstract_declarator = 257, /* direct_abstract_declarator  */
  YYSYMBOL_direct_abstract_declarator_opt = 258, /* direct_abstract_declarator_opt  */
  YYSYMBOL_array_abstract_declarator = 259, /* array_abstract_declarator  */
  YYSYMBOL_function_abstract_declarator = 260, /* function_abstract_declarator  */
  YYSYMBOL_261_9 = 261,                    /* $@9  */
  YYSYMBOL_initializer = 262,              /* initializer  */
  YYSYMBOL_initializer_list = 263,         /* initializer_list  */
  YYSYMBOL_designation_opt = 264,          /* designation_opt  */
  YYSYMBOL_designation = 265,              /* designation  */
  YYSYMBOL_designator_list = 266,          /* designator_list  */
  YYSYMBOL_designator = 267,               /* designator  */
  YYSYMBOL_static_assert_declaration = 268, /* static_assert_declaration  */
  YYSYMBOL_attribute_specifier_sequence = 269, /* attribute_specifier_sequence  */
  YYSYMBOL_attribute_specifier_sequence_opt = 270, /* attribute_specifier_sequence_opt  */
  YYSYMBOL_attribute_specifier = 271,      /* attribute_specifier  */
  YYSYMBOL_attribute_list = 272,           /* attribute_list  */
  YYSYMBOL_attribute = 273,                /* attribute  */
  YYSYMBOL_attribute_opt = 274,            /* attribute_opt  */
  YYSYMBOL_attribute_token = 275,          /* attribute_token  */
  YYSYMBOL_attribute_argument_clause = 276, /* attribute_argument_clause  */
  YYSYMBOL_balanced_token_sequence = 277,  /* balanced_token_sequence  */
  YYSYMBOL_balanced_token = 278,           /* balanced_token  */
  YYSYMBOL_statement = 279,                /* statement  */
  YYSYMBOL_unlabeled_statement = 280,      /* unlabeled_statement  */
  YYSYMBOL_primary_block = 281,            /* primary_block  */
  YYSYMBOL_secondary_block = 282,          /* secondary_block  */
  YYSYMBOL_labeled_statement = 283,        /* labeled_statement  */
  YYSYMBOL_label = 284,                    /* label  */
  YYSYMBOL_285_10 = 285,                   /* @10  */
  YYSYMBOL_start_block = 286,              /* start_block  */
  YYSYMBOL_end_block = 287,                /* end_block  */
  YYSYMBOL_compound_statement = 288,       /* compound_statement  */
  YYSYMBOL_block_item_list = 289,          /* block_item_list  */
  YYSYMBOL_expression_statement = 290,     /* expression_statement  */
  YYSYMBOL_else_statement = 291,           /* else_statement  */
  YYSYMBOL_selection_statement = 292,      /* selection_statement  */
  YYSYMBOL_293_11 = 293,                   /* $@11  */
  YYSYMBOL_294_12 = 294,                   /* @12  */
  YYSYMBOL_iteration_statement = 295,      /* iteration_statement  */
  YYSYMBOL_296_13 = 296,                   /* $@13  */
  YYSYMBOL_297_14 = 297,                   /* $@14  */
  YYSYMBOL_jump_statement = 298,           /* jump_statement  */
  YYSYMBOL_translation_unit = 299,         /* translation_unit  */
  YYSYMBOL_external_declaration = 300,     /* external_declaration  */
  YYSYMBOL_function_definition = 301,      /* function_definition  */
  YYSYMBOL_302_15 = 302,                   /* $@15  */
  YYSYMBOL_303_16 = 303,                   /* $@16  */
  YYSYMBOL_function_body = 304,            /* function_body  */
  YYSYMBOL_file = 305,                     /* file  */
  YYSYMBOL_function_attributes = 306,      /* function_attributes  */
  YYSYMBOL_function_attribute = 307,       /* function_attribute  */
  YYSYMBOL_offsetof_member_designator = 308, /* offsetof_member_designator  */
  YYSYMBOL_309_17 = 309,                   /* $@17  */
  YYSYMBOL_string_literal_val = 310,       /* string_literal_val  */
  YYSYMBOL_Interrupt_storage = 311,        /* Interrupt_storage  */
  YYSYMBOL_sfr_reg_bit = 312,              /* sfr_reg_bit  */
  YYSYMBOL_sfr_attributes = 313,           /* sfr_attributes  */
  YYSYMBOL_opt_stag = 314,                 /* opt_stag  */
  YYSYMBOL_stag = 315,                     /* stag  */
  YYSYMBOL_opt_assign_expr = 316,          /* opt_assign_expr  */
  YYSYMBOL_specifier_qualifier_list = 317, /* specifier_qualifier_list  */
  YYSYMBOL_type_specifier_list_ = 318,     /* type_specifier_list_  */
  YYSYMBOL_identifier_list = 319,          /* identifier_list  */
  YYSYMBOL_type_name = 320,                /* type_name  */
  YYSYMBOL_critical = 321,                 /* critical  */
  YYSYMBOL_critical_statement = 322,       /* critical_statement  */
  YYSYMBOL_statements_and_implicit = 323,  /* statements_and_implicit  */
  YYSYMBOL_declaration_after_statement = 324, /* declaration_after_statement  */
  YYSYMBOL_325_18 = 325,                   /* $@18  */
  YYSYMBOL_implicit_block = 326,           /* implicit_block  */
  YYSYMBOL_declaration_list = 327,         /* declaration_list  */
  YYSYMBOL_kr_declaration = 328,           /* kr_declaration  */
  YYSYMBOL_kr_declaration_list = 329,      /* kr_declaration_list  */
  YYSYMBOL_statement_list = 330,           /* statement_list  */
  YYSYMBOL_while = 331,                    /* while  */
  YYSYMBOL_do = 332,                       /* do  */
  YYSYMBOL_for = 333,                      /* for  */
  YYSYMBOL_asm_string_literal = 334,       /* asm_string_literal  */
  YYSYMBOL_asm_statement = 335,            /* asm_statement  */
  YYSYMBOL_addressmod = 336,               /* addressmod  */
  YYSYMBOL_identifier = 337                /* identifier  */
};
typedef enum yysymbol_kind_t yysymbol_kind_t;




#ifdef short
# undef short
#endif

/* On compilers that do not define __PTRDIFF_MAX__ etc., make sure
   <limits.h> and (if available) <stdint.h> are included
   so that the code can choose integer types of a good width.  */

#ifndef __PTRDIFF_MAX__
# include <limits.h> /* INFRINGES ON USER NAME SPACE */
# if defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stdint.h> /* INFRINGES ON USER NAME SPACE */
#  define YY_STDINT_H
# endif
#endif

/* Narrow types that promote to a signed type and that can represent a
   signed or unsigned integer of at least N bits.  In tables they can
   save space and decrease cache pressure.  Promoting to a signed type
   helps avoid bugs in integer arithmetic.  */

#ifdef __INT_LEAST8_MAX__
typedef __INT_LEAST8_TYPE__ yytype_int8;
#elif defined YY_STDINT_H
typedef int_least8_t yytype_int8;
#else
typedef signed char yytype_int8;
#endif

#ifdef __INT_LEAST16_MAX__
typedef __INT_LEAST16_TYPE__ yytype_int16;
#elif defined YY_STDINT_H
typedef int_least16_t yytype_int16;
#else
typedef short yytype_int16;
#endif

/* Work around bug in HP-UX 11.23, which defines these macros
   incorrectly for preprocessor constants.  This workaround can likely
   be removed in 2023, as HPE has promised support for HP-UX 11.23
   (aka HP-UX 11i v2) only through the end of 2022; see Table 2 of
   <https://h20195.www2.hpe.com/V2/getpdf.aspx/4AA4-7673ENW.pdf>.  */
#ifdef __hpux
# undef UINT_LEAST8_MAX
# undef UINT_LEAST16_MAX
# define UINT_LEAST8_MAX 255
# define UINT_LEAST16_MAX 65535
#endif

#if defined __UINT_LEAST8_MAX__ && __UINT_LEAST8_MAX__ <= __INT_MAX__
typedef __UINT_LEAST8_TYPE__ yytype_uint8;
#elif (!defined __UINT_LEAST8_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST8_MAX <= INT_MAX)
typedef uint_least8_t yytype_uint8;
#elif !defined __UINT_LEAST8_MAX__ && UCHAR_MAX <= INT_MAX
typedef unsigned char yytype_uint8;
#else
typedef short yytype_uint8;
#endif

#if defined __UINT_LEAST16_MAX__ && __UINT_LEAST16_MAX__ <= __INT_MAX__
typedef __UINT_LEAST16_TYPE__ yytype_uint16;
#elif (!defined __UINT_LEAST16_MAX__ && defined YY_STDINT_H \
       && UINT_LEAST16_MAX <= INT_MAX)
typedef uint_least16_t yytype_uint16;
#elif !defined __UINT_LEAST16_MAX__ && USHRT_MAX <= INT_MAX
typedef unsigned short yytype_uint16;
#else
typedef int yytype_uint16;
#endif

#ifndef YYPTRDIFF_T
# if defined __PTRDIFF_TYPE__ && defined __PTRDIFF_MAX__
#  define YYPTRDIFF_T __PTRDIFF_TYPE__
#  define YYPTRDIFF_MAXIMUM __PTRDIFF_MAX__
# elif defined PTRDIFF_MAX
#  ifndef ptrdiff_t
#   include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  endif
#  define YYPTRDIFF_T ptrdiff_t
#  define YYPTRDIFF_MAXIMUM PTRDIFF_MAX
# else
#  define YYPTRDIFF_T long
#  define YYPTRDIFF_MAXIMUM LONG_MAX
# endif
#endif

#ifndef YYSIZE_T
# ifdef __SIZE_TYPE__
#  define YYSIZE_T __SIZE_TYPE__
# elif defined size_t
#  define YYSIZE_T size_t
# elif defined __STDC_VERSION__ && 199901 <= __STDC_VERSION__
#  include <stddef.h> /* INFRINGES ON USER NAME SPACE */
#  define YYSIZE_T size_t
# else
#  define YYSIZE_T unsigned
# endif
#endif

#define YYSIZE_MAXIMUM                                  \
  YY_CAST (YYPTRDIFF_T,                                 \
           (YYPTRDIFF_MAXIMUM < YY_CAST (YYSIZE_T, -1)  \
            ? YYPTRDIFF_MAXIMUM                         \
            : YY_CAST (YYSIZE_T, -1)))

#define YYSIZEOF(X) YY_CAST (YYPTRDIFF_T, sizeof (X))


/* Stored state numbers (used for stacks). */
typedef yytype_int16 yy_state_t;

/* State numbers in computations.  */
typedef int yy_state_fast_t;

#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> /* INFRINGES ON USER NAME SPACE */
#   define YY_(Msgid) dgettext ("bison-runtime", Msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(Msgid) Msgid
# endif
#endif


#ifndef YY_ATTRIBUTE_PURE
# if defined __GNUC__ && 2 < __GNUC__ + (96 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_PURE __attribute__ ((__pure__))
# else
#  define YY_ATTRIBUTE_PURE
# endif
#endif

#ifndef YY_ATTRIBUTE_UNUSED
# if defined __GNUC__ && 2 < __GNUC__ + (7 <= __GNUC_MINOR__)
#  define YY_ATTRIBUTE_UNUSED __attribute__ ((__unused__))
# else
#  define YY_ATTRIBUTE_UNUSED
# endif
#endif

/* Suppress unused-variable warnings by "using" E.  */
#if ! defined lint || defined __GNUC__
# define YY_USE(E) ((void) (E))
#else
# define YY_USE(E) /* empty */
#endif

/* Suppress an incorrect diagnostic about yylval being uninitialized.  */
#if defined __GNUC__ && ! defined __ICC && 406 <= __GNUC__ * 100 + __GNUC_MINOR__
# if __GNUC__ * 100 + __GNUC_MINOR__ < 407
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")
# else
#  define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN                           \
    _Pragma ("GCC diagnostic push")                                     \
    _Pragma ("GCC diagnostic ignored \"-Wuninitialized\"")              \
    _Pragma ("GCC diagnostic ignored \"-Wmaybe-uninitialized\"")
# endif
# define YY_IGNORE_MAYBE_UNINITIALIZED_END      \
    _Pragma ("GCC diagnostic pop")
#else
# define YY_INITIAL_VALUE(Value) Value
#endif
#ifndef YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
# define YY_IGNORE_MAYBE_UNINITIALIZED_END
#endif
#ifndef YY_INITIAL_VALUE
# define YY_INITIAL_VALUE(Value) /* Nothing. */
#endif

#if defined __cplusplus && defined __GNUC__ && ! defined __ICC && 6 <= __GNUC__
# define YY_IGNORE_USELESS_CAST_BEGIN                          \
    _Pragma ("GCC diagnostic push")                            \
    _Pragma ("GCC diagnostic ignored \"-Wuseless-cast\"")
# define YY_IGNORE_USELESS_CAST_END            \
    _Pragma ("GCC diagnostic pop")
#endif
#ifndef YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_BEGIN
# define YY_IGNORE_USELESS_CAST_END
#endif


#define YY_ASSERT(E) ((void) (0 && (E)))

#if !defined yyoverflow

/* The parser invokes alloca or malloc; define the necessary symbols.  */

# ifdef YYSTACK_USE_ALLOCA
#  if YYSTACK_USE_ALLOCA
#   ifdef __GNUC__
#    define YYSTACK_ALLOC __builtin_alloca
#   elif defined __BUILTIN_VA_ARG_INCR
#    include <alloca.h> /* INFRINGES ON USER NAME SPACE */
#   elif defined _AIX
#    define YYSTACK_ALLOC __alloca
#   elif defined _MSC_VER
#    include <malloc.h> /* INFRINGES ON USER NAME SPACE */
#    define alloca _alloca
#   else
#    define YYSTACK_ALLOC alloca
#    if ! defined _ALLOCA_H && ! defined EXIT_SUCCESS
#     include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
      /* Use EXIT_SUCCESS as a witness for stdlib.h.  */
#     ifndef EXIT_SUCCESS
#      define EXIT_SUCCESS 0
#     endif
#    endif
#   endif
#  endif
# endif

# ifdef YYSTACK_ALLOC
   /* Pacify GCC's 'empty if-body' warning.  */
#  define YYSTACK_FREE(Ptr) do { /* empty */; } while (0)
#  ifndef YYSTACK_ALLOC_MAXIMUM
    /* The OS might guarantee only one guard page at the bottom of the stack,
       and a page size can be as small as 4096 bytes.  So we cannot safely
       invoke alloca (N) if N exceeds 4096.  Use a slightly smaller number
       to allow for a few compiler-allocated temporary stack slots.  */
#   define YYSTACK_ALLOC_MAXIMUM 4032 /* reasonable circa 2006 */
#  endif
# else
#  define YYSTACK_ALLOC YYMALLOC
#  define YYSTACK_FREE YYFREE
#  ifndef YYSTACK_ALLOC_MAXIMUM
#   define YYSTACK_ALLOC_MAXIMUM YYSIZE_MAXIMUM
#  endif
#  if (defined __cplusplus && ! defined EXIT_SUCCESS \
       && ! ((defined YYMALLOC || defined malloc) \
             && (defined YYFREE || defined free)))
#   include <stdlib.h> /* INFRINGES ON USER NAME SPACE */
#   ifndef EXIT_SUCCESS
#    define EXIT_SUCCESS 0
#   endif
#  endif
#  ifndef YYMALLOC
#   define YYMALLOC malloc
#   if ! defined malloc && ! defined EXIT_SUCCESS
void *malloc (YYSIZE_T); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
#  ifndef YYFREE
#   define YYFREE free
#   if ! defined free && ! defined EXIT_SUCCESS
void free (void *); /* INFRINGES ON USER NAME SPACE */
#   endif
#  endif
# endif
#endif /* !defined yyoverflow */

#if (! defined yyoverflow \
     && (! defined __cplusplus \
         || (defined YYSTYPE_IS_TRIVIAL && YYSTYPE_IS_TRIVIAL)))

/* A type that is properly aligned for any stack member.  */
union yyalloc
{
  yy_state_t yyss_alloc;
  YYSTYPE yyvs_alloc;
};

/* The size of the maximum gap between one aligned stack and the next.  */
# define YYSTACK_GAP_MAXIMUM (YYSIZEOF (union yyalloc) - 1)

/* The size of an array large to enough to hold all stacks, each with
   N elements.  */
# define YYSTACK_BYTES(N) \
     ((N) * (YYSIZEOF (yy_state_t) + YYSIZEOF (YYSTYPE)) \
      + YYSTACK_GAP_MAXIMUM)

# define YYCOPY_NEEDED 1

/* Relocate STACK from its old location to the new one.  The
   local variables YYSIZE and YYSTACKSIZE give the old and new number of
   elements in the stack, and YYPTR gives the new location of the
   stack.  Advance YYPTR to a properly aligned location for the next
   stack.  */
# define YYSTACK_RELOCATE(Stack_alloc, Stack)                           \
    do                                                                  \
      {                                                                 \
        YYPTRDIFF_T yynewbytes;                                         \
        YYCOPY (&yyptr->Stack_alloc, Stack, yysize);                    \
        Stack = &yyptr->Stack_alloc;                                    \
        yynewbytes = yystacksize * YYSIZEOF (*Stack) + YYSTACK_GAP_MAXIMUM; \
        yyptr += yynewbytes / YYSIZEOF (*yyptr);                        \
      }                                                                 \
    while (0)

#endif

#if defined YYCOPY_NEEDED && YYCOPY_NEEDED
/* Copy COUNT objects from SRC to DST.  The source and destination do
   not overlap.  */
# ifndef YYCOPY
#  if defined __GNUC__ && 1 < __GNUC__
#   define YYCOPY(Dst, Src, Count) \
      __builtin_memcpy (Dst, Src, YY_CAST (YYSIZE_T, (Count)) * sizeof (*(Src)))
#  else
#   define YYCOPY(Dst, Src, Count)              \
      do                                        \
        {                                       \
          YYPTRDIFF_T yyi;                      \
          for (yyi = 0; yyi < (Count); yyi++)   \
            (Dst)[yyi] = (Src)[yyi];            \
        }                                       \
      while (0)
#  endif
# endif
#endif /* !YYCOPY_NEEDED */

/* YYFINAL -- State number of the termination state.  */
#define YYFINAL  154
/* YYLAST -- Last index in YYTABLE.  */
#define YYLAST   3107

/* YYNTOKENS -- Number of terminals.  */
#define YYNTOKENS  185
/* YYNNTS -- Number of nonterminals.  */
#define YYNNTS  153
/* YYNRULES -- Number of rules.  */
#define YYNRULES  410
/* YYNSTATES -- Number of states.  */
#define YYNSTATES  658

/* YYMAXUTOK -- Last valid token kind.  */
#define YYMAXUTOK   415


/* YYTRANSLATE(TOKEN-NUM) -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex, with out-of-bounds checking.  */
#define YYTRANSLATE(YYX)                                \
  (0 <= (YYX) && (YYX) <= YYMAXUTOK                     \
   ? YY_CAST (yysymbol_kind_t, yytranslate[YYX])        \
   : YYSYMBOL_YYUNDEF)

/* YYTRANSLATE[TOKEN-NUM] -- Symbol number corresponding to TOKEN-NUM
   as returned by yylex.  */
static const yytype_uint8 yytranslate[] =
{
       0,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   175,     2,     2,     2,   177,   170,     2,
     161,   162,   171,   172,   163,   173,   167,   176,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,   164,   184,
     178,   183,   179,   182,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,   165,     2,   166,   180,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,   168,   181,   169,   174,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     2,     2,     2,     2,
       2,     2,     2,     2,     2,     2,     1,     2,     3,     4,
       5,     6,     7,     8,     9,    10,    11,    12,    13,    14,
      15,    16,    17,    18,    19,    20,    21,    22,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,    42,    43,    44,
      45,    46,    47,    48,    49,    50,    51,    52,    53,    54,
      55,    56,    57,    58,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    78,    79,    80,    81,    82,    83,    84,
      85,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    96,    97,    98,    99,   100,   101,   102,   103,   104,
     105,   106,   107,   108,   109,   110,   111,   112,   113,   114,
     115,   116,   117,   118,   119,   120,   121,   122,   123,   124,
     125,   126,   127,   128,   129,   130,   131,   132,   133,   134,
     135,   136,   137,   138,   139,   140,   141,   142,   143,   144,
     145,   146,   147,   148,   149,   150,   151,   152,   153,   154,
     155,   156,   157,   158,   159,   160
};

#if YYDEBUG
/* YYRLINE[YYN] -- Source line where rule number YYN was defined.  */
static const yytype_int16 yyrline[] =
{
       0,   174,   174,   175,   176,   177,   178,   179,   183,   184,
     185,   189,   193,   194,   198,   199,   203,   204,   205,   207,
     211,   211,   218,   218,   225,   227,   229,   237,   245,   258,
     259,   263,   264,   265,   266,   275,   276,   277,   278,   279,
     280,   281,   285,   286,   287,   288,   289,   290,   294,   295,
     299,   300,   301,   302,   306,   307,   308,   312,   313,   314,
     318,   319,   320,   321,   322,   326,   327,   328,   332,   333,
     337,   338,   342,   343,   347,   348,   348,   353,   354,   354,
     359,   360,   360,   368,   369,   414,   415,   416,   417,   418,
     419,   420,   421,   422,   423,   424,   428,   429,   429,   433,
     434,   438,   444,   468,   510,   514,   520,   523,   524,   529,
     530,   535,   536,   544,   545,   549,   550,   554,   558,   562,
     566,   570,   575,   579,   583,   590,   595,   600,   605,   610,
     615,   620,   625,   630,   635,   641,   646,   650,   654,   658,
     662,   666,   671,   677,   686,   690,   695,   703,   710,   719,
     722,   727,   732,   743,   757,   756,   867,   897,   898,   902,
     903,   916,   941,   942,   943,   947,   948,   956,   957,   967,
     981,   985,   990,  1014,  1031,  1034,  1042,  1043,  1051,  1079,
    1083,  1087,  1091,  1095,  1099,  1103,  1107,  1111,  1115,  1119,
    1126,  1130,  1137,  1143,  1158,  1159,  1167,  1168,  1169,  1170,
    1174,  1175,  1176,  1180,  1205,  1229,  1254,  1280,  1281,  1313,
    1326,  1325,  1368,  1379,  1380,  1390,  1403,  1423,  1429,  1469,
    1477,  1478,  1486,  1489,  1496,  1497,  1501,  1502,  1510,  1527,
    1537,  1538,  1539,  1546,  1547,  1548,  1552,  1553,  1557,  1564,
    1576,  1577,  1589,  1591,  1590,  1621,  1622,  1630,  1631,  1635,
    1636,  1645,  1646,  1650,  1654,  1655,  1659,  1681,  1685,  1693,
    1706,  1713,  1721,  1724,  1731,  1740,  1744,  1757,  1761,  1769,
    1772,  1779,  1785,  1798,  1799,  1803,  1804,  1808,  1809,  1810,
    1816,  1817,  1822,  1823,  1827,  1831,  1832,  1836,  1840,  1844,
    1851,  1858,  1871,  1876,  1883,  1883,  1893,  1904,  1913,  1914,
    1922,  1923,  1924,  1928,  1929,  1933,  1934,  1939,  1939,  1947,
    1947,  1973,  1973,  1982,  1992,  2018,  2017,  2058,  2068,  2081,
    2091,  2100,  2114,  2115,  2119,  2123,  2144,  2152,  2151,  2171,
    2170,  2201,  2202,  2213,  2216,  2223,  2224,  2228,  2232,  2235,
    2238,  2241,  2244,  2247,  2253,  2256,  2259,  2265,  2271,  2277,
    2280,  2283,  2286,  2289,  2293,  2296,  2299,  2308,  2332,  2350,
    2351,  2351,  2358,  2365,  2392,  2396,  2397,  2412,  2420,  2424,
    2432,  2443,  2454,  2466,  2467,  2478,  2497,  2515,  2528,  2531,
    2533,  2541,  2542,  2550,  2558,  2584,  2595,  2604,  2605,  2614,
    2614,  2625,  2633,  2644,  2660,  2694,  2732,  2744,  2770,  2771,
    2772,  2773,  2776,  2790,  2805,  2830,  2834,  2844,  2857,  2868,
    2886
};
#endif

/** Accessing symbol of state STATE.  */
#define YY_ACCESSING_SYMBOL(State) YY_CAST (yysymbol_kind_t, yystos[State])

#if YYDEBUG || 0
/* The user-facing name of the symbol whose (internal) number is
   YYSYMBOL.  No bounds checking.  */
static const char *yysymbol_name (yysymbol_kind_t yysymbol) YY_ATTRIBUTE_UNUSED;

/* YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
   First, the terminals, then, starting at YYNTOKENS, nonterminals.  */
static const char *const yytname[] =
{
  "\"end of file\"", "error", "\"invalid token\"", "IDENTIFIER",
  "TYPE_NAME", "ADDRSPACE_NAME", "CONSTANT", "SIZEOF", "OFFSETOF",
  "PTR_OP", "INC_OP", "DEC_OP", "LEFT_OP", "RIGHT_OP", "LE_OP", "GE_OP",
  "EQ_OP", "NE_OP", "AND_OP", "OR_OP", "ATTR_START", "TOK_SEP",
  "MUL_ASSIGN", "DIV_ASSIGN", "MOD_ASSIGN", "ADD_ASSIGN", "SUB_ASSIGN",
  "LEFT_ASSIGN", "RIGHT_ASSIGN", "AND_ASSIGN", "XOR_ASSIGN", "OR_ASSIGN",
  "TYPEDEF", "EXTERN", "STATIC", "AUTO", "REGISTER", "CONSTEXPR", "CODE",
  "EEPROM", "INTERRUPT", "SFR", "SFR16", "SFR32", "ADDRESSMOD", "AT",
  "SBIT", "REENTRANT", "USING", "XDATA", "DATA", "IDATA", "PDATA",
  "ELLIPSIS", "CRITICAL", "NONBANKED", "BANKED", "SHADOWREGS", "SD_WPARAM",
  "SD_BOOL", "SD_CHAR", "SD_SHORT", "SD_INT", "SD_LONG", "SIGNED",
  "UNSIGNED", "SD_FLOAT", "DOUBLE", "FIXED16X16", "SD_CONST", "VOLATILE",
  "SD_VOID", "BIT", "COMPLEX", "IMAGINARY", "STRUCT", "UNION", "ENUM",
  "RANGE", "SD_FAR", "CASE", "DEFAULT", "IF", "ELSE", "SWITCH", "WHILE",
  "DO", "FOR", "GOTO", "CONTINUE", "BREAK", "RETURN", "NAKED",
  "JAVANATIVE", "OVERLAY", "TRAP", "STRING_LITERAL", "INLINEASM", "FUNC",
  "IFX", "ADDRESS_OF", "GET_VALUE_AT_ADDRESS", "SET_VALUE_AT_ADDRESS",
  "SPIL", "UNSPIL", "GETABIT", "GETBYTE", "GETWORD", "BITWISEAND",
  "UNARYMINUS", "IPUSH", "IPUSH_VALUE_AT_ADDRESS", "IPOP", "PCALL",
  "ENDFUNCTION", "JUMPTABLE", "RRC", "RLC", "CAST", "CALL", "PARAM",
  "NULLOP", "BLOCK", "LABEL", "RECEIVE", "SEND", "ARRAYINIT",
  "DUMMY_READ_VOLATILE", "ENDCRITICAL", "SWAP", "INLINE", "RESTRICT",
  "SMALLC", "RAISONANCE", "IAR", "COSMIC", "SDCCCALL", "PRESERVES_REGS",
  "Z88DK_FASTCALL", "Z88DK_CALLEE", "Z88DK_SHORTCALL",
  "Z88DK_PARAMS_OFFSET", "ALIGNAS", "ALIGNOF", "ATOMIC", "GENERIC",
  "NORETURN", "STATIC_ASSERT", "THREAD_LOCAL", "TOKEN_FALSE", "TOKEN_TRUE",
  "NULLPTR", "TYPEOF", "TYPEOF_UNQUAL", "SD_BITINT", "DECIMAL32",
  "DECIMAL64", "DECIMAL128", "ASM", "GENERIC_ASSOC_LIST",
  "GENERIC_ASSOCIATION", "'('", "')'", "','", "':'", "'['", "']'", "'.'",
  "'{'", "'}'", "'&'", "'*'", "'+'", "'-'", "'~'", "'!'", "'/'", "'%'",
  "'<'", "'>'", "'^'", "'|'", "'?'", "'='", "';'", "$accept",
  "primary_expression", "predefined_constant", "generic_selection",
  "generic_assoc_list", "generic_association", "postfix_expression", "$@1",
  "$@2", "argument_expr_list", "unary_expression", "unary_operator",
  "cast_expression", "multiplicative_expression", "additive_expression",
  "shift_expression", "relational_expression", "equality_expression",
  "and_expression", "exclusive_or_expression", "inclusive_or_expression",
  "logical_and_expr", "$@3", "logical_or_expr", "$@4", "conditional_expr",
  "$@5", "assignment_expr", "assignment_operator", "expression", "$@6",
  "expression_opt", "constant_expr", "declaration",
  "declaration_specifiers", "declaration_specifiers_",
  "init_declarator_list", "init_declarator", "attribute_declaration",
  "storage_class_specifier", "type_specifier", "typeof_specifier",
  "struct_or_union_specifier", "$@7", "struct_or_union",
  "member_declaration_list", "member_declaration",
  "type_specifier_qualifier", "member_declarator_list",
  "member_declarator", "enum_specifier", "enum_comma_opt",
  "enumerator_list", "enumerator", "type_qualifier", "function_specifier",
  "alignment_specifier", "declarator", "direct_declarator", "declarator2",
  "array_declarator", "declarator2_function_attributes",
  "function_declarator", "$@8", "pointer", "unqualified_pointer",
  "type_qualifier_list", "type_qualifier_list_opt", "parameter_type_list",
  "parameter_list", "parameter_declaration", "abstract_declarator",
  "direct_abstract_declarator", "direct_abstract_declarator_opt",
  "array_abstract_declarator", "function_abstract_declarator", "$@9",
  "initializer", "initializer_list", "designation_opt", "designation",
  "designator_list", "designator", "static_assert_declaration",
  "attribute_specifier_sequence", "attribute_specifier_sequence_opt",
  "attribute_specifier", "attribute_list", "attribute", "attribute_opt",
  "attribute_token", "attribute_argument_clause",
  "balanced_token_sequence", "balanced_token", "statement",
  "unlabeled_statement", "primary_block", "secondary_block",
  "labeled_statement", "label", "@10", "start_block", "end_block",
  "compound_statement", "block_item_list", "expression_statement",
  "else_statement", "selection_statement", "$@11", "@12",
  "iteration_statement", "$@13", "$@14", "jump_statement",
  "translation_unit", "external_declaration", "function_definition",
  "$@15", "$@16", "function_body", "file", "function_attributes",
  "function_attribute", "offsetof_member_designator", "$@17",
  "string_literal_val", "Interrupt_storage", "sfr_reg_bit",
  "sfr_attributes", "opt_stag", "stag", "opt_assign_expr",
  "specifier_qualifier_list", "type_specifier_list_", "identifier_list",
  "type_name", "critical", "critical_statement", "statements_and_implicit",
  "declaration_after_statement", "$@18", "implicit_block",
  "declaration_list", "kr_declaration", "kr_declaration_list",
  "statement_list", "while", "do", "for", "asm_string_literal",
  "asm_statement", "addressmod", "identifier", YY_NULLPTR
};

static const char *
yysymbol_name (yysymbol_kind_t yysymbol)
{
  return yytname[yysymbol];
}
#endif

#define YYPACT_NINF (-513)

#define yypact_value_is_default(Yyn) \
  ((Yyn) == YYPACT_NINF)

#define YYTABLE_NINF (-393)

#define yytable_value_is_error(Yyn) \
  0

/* YYPACT[STATE-NUM] -- Index in YYTABLE of the portion describing
   STATE-NUM.  */
static const yytype_int16 yypact[] =
{
    2228,  -513,  -513,  -513,    63,  -513,  -513,  -513,  -513,  -513,
    -513,  -513,  -513,    25,  -513,  -513,    63,    37,  -513,  -513,
    -513,  -513,  -513,  -513,  -513,  -513,  -513,  -513,  -513,  -513,
    -513,  -513,  -513,  -513,  -513,  -513,  -513,  -513,  -513,  -513,
    -513,    27,  -513,  -513,   -60,  -513,  -513,   -24,  -513,   -21,
      -6,     6,  -513,  -513,  -513,    50,  -513,  -513,    30,  -513,
    -513,  2830,  -513,  -513,  -513,   152,  2830,  -513,  -513,  2830,
    -513,  -513,   -54,    72,     2,   152,  2911,    39,   386,  -513,
       9,  -513,  2228,  -513,  -513,   246,  -513,  -513,  -513,    55,
     -13,  -513,  -513,    97,   240,  -513,    43,  -513,  1967,    63,
     140,  1545,  1967,  1545,  1545,  1967,   124,  -513,   -98,  -513,
     -95,  -513,   152,    63,  -513,  -513,  2458,   542,    11,  -513,
    -513,   151,  -513,   155,  -513,  -513,  -513,  -513,  -513,  -513,
    -513,  -513,  -513,  -513,  -513,  -513,  -513,   157,   159,  -513,
    -513,   161,   172,  2911,  -513,  -513,   -54,   181,  -513,  -513,
    1285,  -513,  -513,  -513,  -513,  -513,    63,   209,    66,  -513,
      63,    63,   164,  -513,  2015,   206,  2079,  2079,  -513,  -513,
    2079,  2079,  2079,   212,   224,  -513,  -513,  -513,  1545,  -513,
    -513,  -513,  -513,  -513,  -513,  -513,  -513,  -513,   105,  -513,
    1967,  -513,   103,   116,   265,    48,   308,   219,   210,   215,
     379,    -2,  -513,   236,  -513,  -513,   242,  -513,   152,    63,
     237,    80,   241,   174,   101,  -513,  -513,   176,   245,   178,
     247,   248,   250,    50,  -513,  1625,  2458,  -513,   244,  -513,
    -513,    50,   633,  -513,  -513,  -513,  2458,   981,   611,  1661,
    -513,  2830,   182,  -513,  1967,  1967,  1967,    63,  1967,  1967,
    -513,  1967,   254,  -513,  -513,  -513,  -513,  -513,  -513,  -513,
      68,  -513,  -513,  -513,   232,  -513,  1545,  -513,  2830,  1545,
    -513,  -513,  -513,  -513,  -513,  2830,  1967,   184,   255,  -513,
    -513,  -513,  1721,  1967,  -513,  -513,  1967,  1967,  1967,  1967,
    1967,  1967,  1967,  1967,  1967,  1967,  1967,  1967,  1967,  1967,
    1967,  1967,  -513,  -513,  -513,  -513,    63,   239,   235,   242,
    -513,  2409,  -111,  -513,   -58,   256,  -513,  -513,  -513,   249,
     323,  -513,  -513,  -513,  -513,  -513,  -513,  -513,  -513,  -513,
    -513,  -513,  1967,  -513,  -513,  -513,  -513,  -513,  -513,  -513,
     258,  -118,  -513,  -513,  -513,   252,   -87,  -513,   259,   262,
    -513,   266,   263,  -513,    30,  1135,   631,  -513,  -513,  -513,
     261,  -513,  1099,  -513,  -513,   808,   920,  -513,   282,  -513,
    -513,   981,  1967,  1967,  -513,   284,    41,   286,   288,  -513,
    -513,  -513,    63,   290,   295,   296,   188,   298,   302,   303,
    1967,  -513,  -513,  -513,   304,   305,   307,   309,   310,  -513,
    1781,    63,  -513,   313,   314,    56,    63,  -513,  -513,  -513,
     103,   103,   116,   116,   265,   265,   265,   265,    48,    48,
     308,   219,   210,  1967,  1967,  1967,  -513,  -513,  1967,  -513,
     301,  -513,   316,   318,   -58,   319,  1818,  -513,   320,  -513,
    1967,  1967,    63,  -513,   -79,  1625,  -513,   -32,  -513,   152,
    -513,  -513,   371,  -513,   -25,  1967,  -513,   315,   322,  -513,
    -513,  -513,    63,   311,   312,  1197,  -513,  -513,  -513,  -513,
    -513,   324,  1099,   325,  -513,  1898,  -513,  -513,  -513,  1099,
    -513,  -513,  -513,  -513,  1018,  2603,  -513,  -513,   321,   331,
    -513,  2278,  -513,    75,  -513,  2652,  -513,  -513,  -513,  -513,
    -513,  1967,  -513,  -513,   327,   332,    63,   332,  -513,  2780,
     118,  -513,  -513,  -513,  1967,  -513,  -513,   215,   379,   189,
    -513,  -513,  -513,  -513,  -513,  2830,  -513,   333,   317,  -513,
     336,  -513,   144,  -513,  -513,  -513,  -513,    28,  -513,  2950,
    -513,   341,  -513,   340,   344,  1967,  1967,   326,  -513,  -513,
    -513,    33,  1967,  -513,   424,  1372,  -513,  -513,  2603,  -513,
    -513,  -513,  -513,   349,  -513,   -26,  -513,   348,   193,  -513,
     350,  -513,   -65,  -513,  1967,   353,  -513,  -513,  -513,  -513,
    1625,  -513,  -513,  -513,    54,  2950,   334,  -513,  -513,   195,
     197,  -513,  -513,   201,   355,   335,  1967,  -513,  -513,  1967,
    -513,  1967,  -513,  2780,  1967,   150,  -513,  -513,  -513,  -513,
    1967,    36,  -513,   357,  -513,  -513,  -513,  -513,  -513,  1967,
    1967,   338,    93,    63,  -513,  -513,  -513,  -513,  -513,    54,
    -513,  1967,  1099,  1099,  1099,   216,   339,  1967,  -513,  -513,
    -513,  -513,   441,  -513,  -513,   342,  1967,   363,  1099,  -513,
    -513,   365,  -513,  -513,  1099,  1099,  -513,  -513
};

/* YYDEFACT[STATE-NUM] -- Default reduction number in state STATE-NUM.
   Performed when YYTABLE does not specify something else to do.  Zero
   means the default is an error.  */
static const yytype_int16 yydefact[] =
{
     333,   410,   143,   183,   269,   118,   119,   120,   122,   123,
     124,   185,   186,   369,   371,   372,     0,     0,   367,   184,
     187,   188,   189,   135,   126,   127,   128,   129,   132,   133,
     130,   131,   145,   179,   181,   125,   146,   136,   137,   157,
     158,     0,   190,   180,     0,   182,   191,     0,   121,     0,
       0,     0,   138,   139,   140,     0,   219,   325,     0,   106,
     105,   107,   162,   144,   141,   262,   109,   142,   163,   111,
     164,   327,   194,     0,   262,   262,   207,     0,   213,   104,
       0,   261,   334,   322,   324,     0,   149,   368,   326,   262,
       0,   270,   265,   267,   271,   370,     0,   147,     0,     0,
     173,     0,     0,     0,     0,     0,     0,   102,     0,   113,
     329,   108,   263,   374,   110,   112,     0,   222,   210,   198,
     199,   365,   338,     0,   339,   343,   346,   344,   345,   340,
     341,   342,   348,   349,   350,   351,   352,     0,     0,   354,
     355,     0,     0,   208,   335,   347,   195,     0,   220,   217,
     215,   117,   260,   323,     1,   196,   269,     0,     0,   268,
       0,     0,     0,     3,     0,     0,     0,     0,   364,   363,
       0,     0,     0,     0,     0,     8,     9,    10,     0,    42,
      43,    44,    45,    46,    47,    16,     7,     6,    31,    48,
       0,    50,    54,    57,    60,    65,    68,    70,    72,    74,
      77,    80,   101,     0,     4,     2,   174,   176,   262,     0,
       0,   383,     0,     0,    48,    83,    96,     0,     0,     0,
       0,     0,   197,     0,   103,     0,     0,   154,   156,   375,
     296,     0,   262,   331,   328,   396,     0,   222,   223,     0,
     209,     0,     0,   381,     0,     0,     0,     0,     0,     0,
     336,     0,     0,   221,   218,   266,   264,   279,   278,   273,
       0,   275,   277,   272,     0,   408,     0,    35,     0,     0,
      32,    33,    40,    39,    41,     0,     0,     0,     0,    22,
      24,    25,     0,     0,    20,    34,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    75,    78,    81,   148,   175,     0,   377,   174,
     193,   236,   230,   384,   231,     0,   234,   235,   192,     0,
       0,    86,    87,    88,    89,    90,    91,    92,    93,    94,
      95,    85,     0,   150,    97,   151,   152,   153,   134,   114,
     115,   251,   245,   116,   330,     0,     0,   385,     0,     0,
     297,   100,     0,   393,     0,   263,     0,   398,   399,   298,
       0,   282,   262,   285,   300,   262,   389,   286,     2,   332,
     397,   223,     0,     0,   203,     0,   383,     0,   224,   226,
     229,   212,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   274,   276,   409,     0,     0,     0,     0,     0,     5,
       0,     0,    18,     0,    29,     0,     0,    51,    52,    53,
      55,    56,    58,    59,    63,    64,    61,    62,    66,    67,
      69,    71,    73,     0,     0,     0,   177,   171,     0,   178,
       0,   240,     0,     0,   232,   243,     0,   259,     0,    84,
       0,     0,     0,   246,     0,     0,   252,     0,   254,   262,
     395,   407,     0,   303,     0,     0,   294,     0,     0,   402,
     403,   404,     0,     0,     0,     0,   283,   287,   288,   289,
     284,     0,   262,     0,   299,   263,   386,   281,   280,   262,
     394,   302,   400,   401,   262,     0,   388,   292,     0,     0,
     204,   236,   228,   230,   211,     0,   382,   366,   337,   353,
     358,     0,   356,   214,     0,    36,     0,     0,    37,     0,
     251,    49,    23,    19,     0,    17,    21,    76,    79,     0,
     376,   172,   242,   233,   241,     0,   238,     0,     0,    98,
       0,   257,   251,   247,   249,   253,   255,   262,   159,     0,
     405,     0,   304,     0,     0,     0,     0,     0,   318,   319,
     320,     0,     0,   290,     0,    99,   291,   391,   390,   205,
     206,   225,   227,     0,   216,     0,   359,     0,     0,    12,
       0,    28,     0,    30,     0,     0,   239,   258,   256,   248,
       0,   155,   160,   379,   170,   378,     0,   293,   295,     0,
       0,   317,   321,     0,     0,     0,    99,   357,    38,     0,
     360,     0,    11,     0,     0,   251,    26,    82,   244,   250,
       0,     0,   165,   167,   380,   406,   307,   309,   311,     0,
      99,     0,     0,     0,    15,    13,    14,    27,   168,   170,
     161,     0,   262,   262,   262,     0,     0,    99,   362,   361,
     166,   169,   306,   310,   312,     0,    99,     0,   262,   308,
     313,     0,   315,   305,   262,   262,   314,   316
};

/* YYPGOTO[NTERM-NUM].  */
static const yytype_int16 yypgoto[] =
{
    -513,  -513,  -513,  -513,  -513,   -74,  -513,  -513,  -513,    17,
       4,  -513,  -166,    91,    92,    34,    95,   233,   238,   234,
     113,   115,  -513,  -513,  -513,   -67,  -513,    -1,  -513,   -12,
    -513,  -503,   -41,  -223,     3,   207,   306,   328,  -513,  -513,
    -513,  -513,  -513,  -513,  -513,  -513,     7,  -512,  -513,   -88,
    -513,   243,   337,   253,   -50,  -513,  -513,     1,   -71,  -513,
    -513,  -513,  -513,  -513,   -42,  -513,   -85,   329,  -218,  -513,
      53,  -301,  -304,  -513,  -513,  -513,  -513,  -426,    35,  -495,
    -513,  -513,   102,  -513,     5,   -63,   -73,  -513,  -513,   394,
    -513,  -513,  -513,   293,  -322,  -219,  -513,  -391,  -513,  -217,
    -513,  -513,   194,  -191,  -513,  -513,  -513,  -513,  -513,  -513,
    -513,  -513,  -513,  -513,  -513,   473,  -513,  -513,  -513,   330,
    -513,  -513,   419,  -513,  -513,  -513,  -513,  -513,  -513,  -513,
    -513,  -513,  -513,  -513,   351,   -83,  -513,  -513,  -340,  -513,
    -513,  -513,    82,   346,  -513,  -513,  -513,  -513,  -513,  -513,
    -513,  -513,     0
};

/* YYDEFGOTO[NTERM-NUM].  */
static const yytype_int16 yydefgoto[] =
{
       0,   185,   186,   187,   568,   569,   188,   406,   401,   403,
     214,   190,   191,   192,   193,   194,   195,   196,   197,   198,
     199,   200,   423,   201,   424,   215,   425,   216,   332,   351,
     440,   352,   203,    57,   211,    59,   108,   109,    60,    61,
      62,    63,    64,   345,    65,   537,   538,    66,   611,   612,
      67,   307,   206,   207,    68,    69,    70,   340,    72,    73,
      74,    75,    76,   241,    77,    78,   150,   239,   432,   378,
     379,   313,   314,   315,   316,   317,   525,   343,   444,   445,
     446,   447,   448,    79,   475,   356,    81,    90,    91,    92,
      93,   159,   260,   261,   553,   477,   466,   554,   478,   479,
     544,   232,   359,   233,   360,   361,   649,   468,   632,   633,
     469,   634,   655,   470,    82,    83,    84,   116,   226,   234,
      85,   143,   144,   565,   623,   204,   145,    86,    87,   227,
     228,   429,   584,   585,   242,   380,   362,   363,   364,   484,
     485,   486,   365,   235,   236,   366,   471,   472,   473,   541,
     367,    88,   205
};

/* YYTABLE[YYPACT[STATE-NUM]] -- What to do in state STATE-NUM.  If
   positive, shift that token.  If negative, reduce the rule whose
   number is the opposite.  If YYTABLE_NINF, syntax error.  */
static const yytype_int16 yytable[] =
{
      89,    71,   113,    58,    94,    80,   146,   152,   434,   353,
     433,   119,   120,   357,     1,   358,    96,   303,   212,   534,
     218,   220,     4,   377,   285,   481,   155,   583,   148,     4,
       1,   202,   238,     1,   202,   202,   149,   580,   202,   152,
     476,   100,     1,    97,     1,   369,     1,   441,     4,   442,
     311,   443,   595,     1,  -236,    89,   106,     1,    89,   110,
     210,   213,   293,   294,   221,   223,     1,   148,  -115,     1,
     112,     1,   257,   614,   257,     4,   223,    89,     1,   112,
     112,    95,    89,    71,   532,    58,   224,    80,   225,  -115,
     533,   217,   219,   621,   112,   278,   162,   450,   605,   208,
     253,   101,   189,   435,   606,   189,   189,  -237,   254,   189,
     580,   117,   161,   229,   279,   280,   281,   636,   243,   231,
     407,   408,   409,   321,   322,   323,   324,   325,   326,   327,
     328,   329,   330,   441,   647,   442,   598,   102,   334,   599,
     103,   600,   480,   651,   557,   308,   357,   482,   358,   483,
     156,   535,   371,   157,   609,   104,    94,   556,   262,   542,
     263,   264,   258,  -202,   258,   467,   277,   105,   267,   312,
     270,   271,     4,   240,   272,   273,   274,   202,   202,   202,
     304,   202,   202,   394,   202,   395,   396,   148,   253,   434,
     433,    55,   397,   151,   189,    99,   334,   581,    98,   629,
      55,    56,   491,   383,   384,   385,  -236,   387,   388,   208,
     389,    55,    56,   112,   107,    55,  -200,   592,   610,   334,
     630,    56,   515,    89,   342,    56,   295,   296,   259,   231,
     391,    89,   368,   118,   511,   354,   491,   355,   375,   231,
    -236,   311,   643,   644,   376,  -236,   154,   243,   189,   189,
     189,    56,   189,   189,   277,   189,   334,   277,   158,   638,
     262,   160,   353,   656,   657,   357,   282,   358,   111,   312,
     283,   405,   284,   114,   286,   398,   115,   291,   292,   287,
     288,   404,   152,   441,   331,   442,   222,   571,   289,   290,
     189,   189,   189,   189,   189,   189,   189,   189,   189,   189,
     189,   189,   189,   189,   189,   189,   208,   575,   209,   441,
     642,   442,   244,   579,   376,   441,   245,   442,   246,   627,
     247,   253,   248,   202,   297,   298,   653,   414,   415,   416,
     417,   439,   596,   249,   493,   480,   319,   320,   333,   334,
     336,   334,   251,   454,   381,   382,   399,   334,   265,   504,
     500,   382,   334,   574,    89,   602,   603,   616,   334,   617,
     334,   202,   368,   618,   334,   368,   368,   268,   354,   202,
     355,   488,   489,   275,   202,   256,    89,   492,   645,   334,
     410,   411,   496,   412,   413,   276,   539,   520,   202,   299,
     300,     3,   418,   419,   189,   527,   301,   302,   305,   310,
     530,   512,   152,   318,   189,   306,   516,   335,   427,   337,
     338,  -201,  -373,   519,   543,   390,   393,   400,   428,   438,
     449,   436,   146,   452,    11,    12,   570,   189,   189,   334,
     350,   147,   189,   437,   202,    19,    20,    21,    22,   529,
     189,   225,   531,   451,   342,   189,   487,   453,   494,   493,
     490,   495,   497,   551,   112,    33,    34,   498,   499,   189,
     563,   501,   547,   454,   502,   503,   505,   540,   506,   507,
     521,   508,   368,   509,   539,   513,   545,   514,   522,   368,
     523,   524,   528,   546,   368,   552,   555,   559,   354,   564,
      80,    89,   106,    89,   376,   548,   549,   560,   376,   576,
     510,   577,   578,   586,   587,   189,   566,   607,   588,   594,
     591,   597,   601,   404,   604,   608,   619,    43,   615,   620,
     570,   631,   637,   646,   648,   652,   650,   654,   376,   625,
      45,   573,   420,   589,   590,   422,   517,   346,   421,   518,
     593,   640,   112,   202,   582,   572,   309,     3,   562,   536,
     255,   339,   430,   392,   474,   153,   344,    56,   354,   426,
      80,   354,   250,    80,   202,     0,   372,   558,     0,   628,
       0,     0,     0,     0,     0,     0,   237,     0,   189,   342,
      11,    12,   370,     0,    89,   613,     0,   622,     0,     0,
     641,    19,    20,    21,    22,     0,     0,     0,   386,     0,
     624,     0,     0,   626,     0,     0,     0,   635,     0,     0,
       0,    33,    34,     0,   189,     0,     3,     0,     0,     0,
       0,     0,     0,   639,     0,     0,     0,     0,     0,    89,
     613,     0,   368,   368,   368,   189,     1,     2,     3,   163,
     164,   165,     0,   166,   167,   373,     0,     0,   368,    11,
      12,     0,     0,     4,   368,   368,     0,     0,     0,     0,
      19,    20,    21,    22,     0,     5,     6,     7,     8,     9,
      10,    11,    12,    43,    13,    14,    15,     0,    17,    18,
      33,    34,    19,    20,    21,    22,    45,   347,     0,     0,
       0,     0,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,   455,   456,   457,     0,   458,   459,   460,   461,   462,
     463,   464,   465,     0,     0,     0,     0,     0,     0,   168,
     348,   169,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    43,     0,     0,     0,     0,     0,     0,   170,
     171,     0,     0,     0,     0,    45,     0,     0,     0,     0,
       0,     0,   172,    42,    43,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    44,   173,    45,   174,    46,
      47,    48,   175,   176,   177,    49,    50,    51,    52,    53,
      54,   349,     0,     0,   178,     0,     0,     0,     0,   230,
       0,     0,   350,   179,   180,   181,   182,   183,   184,     0,
       0,     1,     2,     3,   163,   164,   165,   -99,   166,   167,
       0,     0,     0,     0,     0,     0,     0,     0,     4,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       5,     6,     7,     8,     9,    10,    11,    12,     0,    13,
      14,    15,     0,    17,    18,     0,     0,    19,    20,    21,
      22,     0,   347,     0,     0,     0,     0,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   168,   348,   169,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     1,   170,   171,   163,   164,   165,     0,
     166,   167,     0,     0,     0,     0,     0,   172,    42,    43,
       4,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      44,   173,    45,   174,    46,    47,    48,   175,   176,   177,
      49,    50,    51,    52,    53,    54,   349,     0,     0,   178,
       0,     0,     0,     0,   347,     0,     0,  -301,   179,   180,
     181,   182,   183,   184,     0,     0,     3,     0,     0,     0,
       0,     0,   -99,     0,     0,     0,     0,     0,     0,     0,
    -262,  -262,  -262,     0,  -262,  -262,  -262,  -262,  -262,  -262,
    -262,  -262,     0,     0,     0,     0,   168,   348,   169,    11,
      12,     1,     0,     0,   163,   164,   165,     0,   166,   167,
      19,    20,    21,    22,     0,     0,   170,   171,     4,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   172,
      33,    34,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   173,     0,   174,     0,     0,     0,   175,
     176,   177,   347,     0,     0,     0,     0,     0,   349,     0,
       0,   178,     0,     0,     0,     0,     0,     0,  -262,  -387,
     179,   180,   181,   182,   183,   184,     0,     0,     0,     0,
       0,     0,     1,     0,   -99,   163,   164,   165,     0,   166,
     167,     0,    43,     0,   168,   348,   169,     0,     0,     4,
       0,     0,     0,     0,     0,    45,     0,     0,     0,     0,
       0,     0,     0,     0,   170,   171,     0,     0,     1,     0,
       0,   163,   164,   165,     0,   166,   167,   172,     0,     0,
       0,     0,     0,   347,     0,     4,     0,     0,     0,     0,
       0,   173,     0,   174,     0,     0,     0,   175,   176,   177,
       0,     0,     0,     0,     0,     0,   349,     0,     0,   178,
       0,     0,     0,     0,     0,     0,     0,  -392,   179,   180,
     181,   182,   183,   184,     0,   168,   348,   169,     0,     0,
       1,     0,   -99,   163,   164,   165,     0,   166,   167,     0,
       0,     0,     0,     0,     0,   170,   171,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   172,     0,
       0,   168,     0,   169,     0,     0,     0,     0,     0,     0,
       0,     0,   173,     0,   174,     0,     0,     0,   175,   176,
     177,   170,   171,     0,     0,     0,     0,   349,     0,     0,
     178,     0,     0,     0,   172,     0,     0,     0,     0,   179,
     180,   181,   182,   183,   184,     0,     0,     0,   173,     0,
     174,     0,     0,   -99,   175,   176,   177,     0,     0,     0,
       3,     0,     0,   168,     0,   169,   178,     0,     0,     0,
       0,     0,     0,     0,     0,   179,   180,   181,   182,   183,
     184,     0,     0,   170,   171,     0,     0,     0,     0,   151,
       0,     0,     0,    11,    12,     0,   172,     0,     0,     0,
     252,     0,     0,     0,    19,    20,    21,    22,     0,     0,
     173,     0,   174,     0,     0,     0,   175,   176,   177,     0,
       0,     0,     0,     0,    33,    34,     0,     0,   178,     0,
       0,     0,     0,     0,     0,     0,     0,   179,   180,   181,
     182,   183,   184,     0,     0,     1,     2,     3,   163,   164,
     165,   550,   166,   167,     0,     0,     0,     0,     0,     0,
       0,     0,     4,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     5,     6,     7,     8,     9,    10,
      11,    12,     0,    13,    14,    15,    43,    17,    18,     0,
       0,    19,    20,    21,    22,     0,     0,     0,     0,    45,
       0,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
       0,     0,     0,     0,     0,     0,    56,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   168,     0,
     169,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   170,   171,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   172,    42,    43,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,    44,   173,    45,   174,    46,    47,
      48,   175,   176,   177,    49,    50,    51,    52,    53,    54,
       0,     0,     0,   178,     0,     0,     0,     0,     0,     0,
       0,     0,   179,   180,   181,   182,   183,   184,     1,     2,
       3,   163,   164,   165,     0,   166,   167,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     5,     6,     7,
       8,     9,    10,    11,    12,     0,    13,    14,    15,     0,
      17,    18,     0,     0,    19,    20,    21,    22,     0,     0,
       0,     0,     0,     0,    23,    24,    25,    26,    27,    28,
      29,    30,    31,    32,    33,    34,    35,    36,    37,    38,
      39,    40,    41,     0,     0,     0,     0,     0,     1,     0,
       0,   163,   164,   165,     0,   166,   167,     0,     0,     0,
       0,   168,     0,   169,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   170,   171,     0,     1,     0,     0,   163,   164,   165,
       0,   166,   167,     0,   172,    42,    43,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    44,   173,    45,
     174,    46,     0,    48,   175,   176,   177,    49,    50,    51,
      52,    53,    54,     0,     0,     0,   178,     0,     0,     0,
       0,     0,     0,     0,     0,   179,   180,   181,   182,   183,
     184,   168,     0,   169,     1,     0,     0,   163,   164,   165,
       0,   166,   167,     0,     0,     0,     0,     0,     0,     0,
       0,   170,   171,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   172,     0,     0,   168,     0,   169,
       0,     0,     0,     0,     0,     0,     0,     0,   173,     0,
     174,     0,     0,     0,   175,   176,   177,   170,   171,     0,
       0,     0,     0,     0,     1,     0,   178,   163,   164,   165,
     172,   166,   167,   341,     0,   179,   180,   181,   182,   183,
     184,     0,     0,     0,   173,     0,   174,     0,     0,     0,
     175,   176,   177,     0,     0,     0,     0,   168,     0,   169,
       0,     1,   178,     0,   163,   164,   165,   374,   166,   167,
       0,   179,   180,   181,   182,   183,   184,   170,   171,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     172,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   173,     0,   174,     0,     0,     0,
     175,   176,   177,     0,     0,     0,     0,   168,     0,   169,
       0,     0,   178,   402,     0,     0,     0,     0,     0,     0,
       0,   179,   180,   181,   182,   183,   184,   170,   171,     0,
       0,     1,     0,     0,   163,   164,   165,     0,   166,   167,
     172,     0,     0,     0,   168,     0,   169,     0,     4,     0,
       0,     0,     0,     0,   173,     0,   174,     0,     0,     0,
     175,   176,   177,     0,   170,   171,     0,     0,     0,     0,
       0,     0,   178,     0,     0,     0,     0,   172,     0,   510,
       0,   179,   180,   181,   182,   183,   184,     0,     0,     0,
       0,   173,     0,   174,     0,     0,     0,   175,   176,   177,
       1,     0,     0,   163,   164,   165,     0,   166,   167,   178,
       0,     0,     0,     0,   526,     0,     0,     0,   179,   180,
     181,   182,   183,   184,   168,     0,   169,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,   170,   171,     0,     0,     1,     0,
       0,   163,   164,   165,     0,   166,   167,   172,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   173,     0,   174,     0,     0,     0,   175,   176,   177,
       0,     0,     0,     0,     0,     0,     0,     0,     0,   178,
       0,     0,     0,   168,     0,   169,     0,     0,   179,   180,
     181,   182,   183,   184,     0,     0,     0,     0,     0,     0,
       0,     0,     1,   170,   171,   163,   164,   165,     0,   166,
     167,     0,     0,     0,     0,     0,   172,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
     173,   168,   174,   169,     0,     0,   175,   176,   177,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   178,     0,
       0,   170,   171,     0,     0,     0,     0,   179,   180,   181,
     182,   183,   184,     0,   172,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   173,     0,
     174,     0,     0,     0,   175,   176,   177,     0,     0,     0,
       0,     0,     0,     0,     0,   168,   266,   169,     0,     0,
       0,     0,     0,     0,     0,   179,   180,   181,   182,   183,
     184,     0,     0,     0,     0,   170,   171,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,   172,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,   173,     0,   174,     0,     0,     0,   175,   176,
     177,     1,     2,     3,     0,     0,     0,     0,     0,     0,
     269,     0,     0,     0,     0,     0,     0,     0,     4,   179,
     180,   181,   182,   183,   184,     0,     0,     0,     0,     0,
       5,     6,     7,     8,     9,    10,    11,    12,     0,    13,
      14,    15,    16,    17,    18,     0,     0,    19,    20,    21,
      22,     1,     2,     3,     0,     0,     0,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,     0,     0,     0,     0,
       5,     6,     7,     8,     9,    10,    11,    12,     0,    13,
      14,    15,     0,    17,    18,     0,     0,    19,    20,    21,
      22,     0,     0,     0,     0,     0,     0,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,     0,     0,    42,    43,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      44,     0,    45,     0,    46,    47,    48,     0,     0,     0,
      49,    50,    51,    52,    53,    54,     0,     0,     0,    55,
       0,     0,     0,     0,     0,     0,     0,     0,     0,    56,
       0,     0,     0,     0,     0,     0,     0,     0,    42,    43,
       0,     0,     0,     2,     3,     0,     0,     0,     0,     0,
      44,     0,    45,     0,    46,     0,    48,     0,     0,     0,
      49,    50,    51,    52,    53,    54,     0,     0,     0,   491,
     431,     5,     6,     7,     8,     9,    10,    11,    12,    56,
      13,    14,    15,     0,    17,    18,     0,     0,    19,    20,
      21,    22,     2,     3,     0,     0,     0,     0,    23,    24,
      25,    26,    27,    28,    29,    30,    31,    32,    33,    34,
      35,    36,    37,    38,    39,    40,    41,     0,     0,     0,
       5,     6,     7,     8,     9,    10,    11,    12,     0,    13,
      14,    15,     0,    17,    18,     0,     0,    19,    20,    21,
      22,     0,     0,     0,     0,     0,     0,    23,    24,    25,
      26,    27,    28,    29,    30,    31,    32,    33,    34,    35,
      36,    37,    38,    39,    40,    41,     0,     0,     0,    42,
      43,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    44,     0,    45,     0,    46,     0,    48,     0,     0,
       0,    49,    50,    51,    52,    53,    54,     0,     0,     0,
     311,   431,     0,     0,     0,     0,     0,     0,     0,     0,
      56,     0,     0,     0,     0,     0,     0,     0,    42,    43,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
      44,     0,    45,     0,    46,     0,    48,     2,     3,     0,
      49,    50,    51,    52,    53,    54,     0,     0,     0,     0,
       0,     0,     0,     4,     0,     0,   230,     0,     0,     0,
       0,     0,     0,     0,     0,     5,     6,     7,     8,     9,
      10,    11,    12,     0,    13,    14,    15,     0,    17,    18,
       0,     0,    19,    20,    21,    22,     2,     3,     0,     0,
       0,     0,    23,    24,    25,    26,    27,    28,    29,    30,
      31,    32,    33,    34,    35,    36,    37,    38,    39,    40,
      41,     0,     0,     0,     5,     6,     7,     8,     9,    10,
      11,    12,     0,    13,    14,    15,     0,    17,    18,     0,
       0,    19,    20,    21,    22,   561,     0,     0,     0,     0,
       0,    23,    24,    25,    26,    27,    28,    29,    30,    31,
      32,    33,    34,    35,    36,    37,    38,    39,    40,    41,
       0,     0,     0,    42,    43,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,    44,     0,    45,     0,    46,
      47,    48,     0,     0,     0,    49,    50,    51,    52,    53,
      54,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    42,    43,     2,     3,     0,     0,     0,     0,
       0,     0,     0,     0,    44,     0,    45,     0,    46,     0,
      48,     0,     0,     0,    49,    50,    51,    52,    53,    54,
       0,     0,     5,     6,     7,     8,     9,    10,    11,    12,
       0,    13,    14,    15,     0,    17,    18,     0,     0,    19,
      20,    21,    22,     0,     2,     3,     0,     0,     0,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,     0,     0,
       0,   567,     5,     6,     7,     8,     9,    10,    11,    12,
       0,    13,    14,    15,     0,    17,    18,     0,     0,    19,
      20,    21,    22,     0,     0,     0,     0,     0,     0,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,     0,     0,
      42,    43,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    44,     0,    45,     0,    46,     0,    48,     0,
       0,     0,    49,    50,    51,    52,    53,    54,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,   121,     0,     0,     2,     3,     0,     0,   122,   123,
      42,    43,     0,     0,     0,   124,   125,   126,   127,   128,
       0,     0,    44,     0,    45,     0,    46,     0,    48,     0,
       0,     0,    49,    50,    51,    52,    53,    54,    11,    12,
       0,    13,    14,    15,     0,    17,    18,     0,     0,    19,
      20,    21,    22,   129,   130,   131,   132,     0,     0,    23,
      24,    25,    26,    27,    28,    29,    30,    31,    32,    33,
      34,    35,    36,    37,    38,    39,    40,    41,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,   133,   134,   135,   136,   137,   138,   139,
     140,   141,   142,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,     0,     0,     0,
       0,    43,     0,     0,     0,     0,     0,     0,     0,     0,
       0,     0,    44,     0,    45,     0,     0,     0,     0,     0,
       0,     0,    49,    50,    51,    52,    53,    54
};

static const yytype_int16 yycheck[] =
{
       0,     0,    65,     0,     4,     0,    77,    80,   312,   232,
     311,    74,    75,   232,     3,   232,    16,    19,   101,   445,
     103,   104,    20,   241,   190,   365,    89,   539,    78,    20,
       3,    98,   117,     3,   101,   102,    78,   532,   105,   112,
     362,    41,     3,     6,     3,   236,     3,   165,    20,   167,
     161,   169,   555,     3,   165,    55,    55,     3,    58,    58,
     101,   102,    14,    15,   105,   163,     3,   117,   163,     3,
      65,     3,     6,   585,     6,    20,   163,    77,     3,    74,
      75,    56,    82,    82,   163,    82,   184,    82,   183,   184,
     169,   103,   104,   596,    89,   178,    96,   184,   163,    99,
     150,   161,    98,   161,   169,   101,   102,   165,   150,   105,
     605,   165,    69,   113,     9,    10,    11,   620,   118,   116,
     286,   287,   288,    22,    23,    24,    25,    26,    27,    28,
      29,    30,    31,   165,   637,   167,   162,   161,   163,   165,
     161,   167,   365,   646,   484,   208,   365,   366,   365,   366,
     163,   183,   237,   166,   580,   161,   156,   479,   158,   184,
     160,   161,    96,   161,    96,   356,   178,   161,   164,   211,
     166,   167,    20,   162,   170,   171,   172,   244,   245,   246,
     182,   248,   249,   266,   251,   268,   269,   237,   238,   493,
     491,   161,   275,   184,   190,   168,   163,   169,   161,   163,
     161,   171,   161,   244,   245,   246,   165,   248,   249,   209,
     251,   161,   171,   208,   184,   161,   161,   184,   164,   163,
     184,   171,   166,   223,   225,   171,   178,   179,   162,   226,
     162,   231,   232,   161,   400,   232,   161,   232,   239,   236,
     165,   161,   633,   634,   241,   165,     0,   247,   244,   245,
     246,   171,   248,   249,   266,   251,   163,   269,   161,   166,
     260,    21,   485,   654,   655,   484,   161,   484,    61,   311,
     165,   283,   167,    66,   171,   276,    69,    12,    13,   176,
     177,   282,   355,   165,   183,   167,   162,   169,   172,   173,
     286,   287,   288,   289,   290,   291,   292,   293,   294,   295,
     296,   297,   298,   299,   300,   301,   306,   525,   168,   165,
     632,   167,   161,   169,   311,   165,   161,   167,   161,   169,
     161,   371,   161,   390,    16,    17,   648,   293,   294,   295,
     296,   332,   555,   161,   376,   558,   162,   163,   162,   163,
     162,   163,   161,   355,   162,   163,   162,   163,   184,   390,
     162,   163,   163,   164,   354,   162,   163,   162,   163,   162,
     163,   428,   362,   162,   163,   365,   366,   161,   365,   436,
     365,   372,   373,   161,   441,   166,   376,   376,   162,   163,
     289,   290,   382,   291,   292,   161,   449,   428,   455,   170,
     180,     5,   297,   298,   390,   436,   181,    18,   162,   162,
     441,   401,   475,   162,   400,   163,   406,   162,   169,   162,
     162,   161,   168,   425,   455,   161,   184,   162,   183,    96,
     168,   165,   493,   161,    38,    39,   509,   423,   424,   163,
     169,    45,   428,   184,   501,    49,    50,    51,    52,   440,
     436,   183,   442,   184,   445,   441,   164,   184,   162,   491,
     166,   163,   162,   465,   449,    69,    70,   162,   162,   455,
     501,   163,   462,   475,   162,   162,   162,    96,   163,   162,
     169,   162,   472,   163,   537,   162,   161,   163,   162,   479,
     162,   162,   162,   161,   484,   161,   161,   166,   485,   162,
     485,   491,   491,   493,   491,   184,   184,   166,   495,   166,
     168,   184,   166,   162,   164,   501,   506,   574,   164,    85,
     184,   162,   164,   514,   164,   162,   161,   131,   184,   184,
     603,   164,   184,   184,    83,   162,   184,   162,   525,   603,
     144,   514,   299,   545,   546,   301,   423,   231,   300,   424,
     552,   629,   537,   610,   537,   510,   209,     5,   495,   447,
     156,   223,   309,   260,   360,    82,   226,   171,   555,   306,
     555,   558,   143,   558,   631,    -1,   237,   485,    -1,   610,
      -1,    -1,    -1,    -1,    -1,    -1,    34,    -1,   574,   580,
      38,    39,   236,    -1,   584,   584,    -1,   599,    -1,    -1,
     631,    49,    50,    51,    52,    -1,    -1,    -1,   247,    -1,
     601,    -1,    -1,   604,    -1,    -1,    -1,   619,    -1,    -1,
      -1,    69,    70,    -1,   610,    -1,     5,    -1,    -1,    -1,
      -1,    -1,    -1,   623,    -1,    -1,    -1,    -1,    -1,   629,
     629,    -1,   632,   633,   634,   631,     3,     4,     5,     6,
       7,     8,    -1,    10,    11,    34,    -1,    -1,   648,    38,
      39,    -1,    -1,    20,   654,   655,    -1,    -1,    -1,    -1,
      49,    50,    51,    52,    -1,    32,    33,    34,    35,    36,
      37,    38,    39,   131,    41,    42,    43,    -1,    45,    46,
      69,    70,    49,    50,    51,    52,   144,    54,    -1,    -1,
      -1,    -1,    59,    60,    61,    62,    63,    64,    65,    66,
      67,    68,    69,    70,    71,    72,    73,    74,    75,    76,
      77,    80,    81,    82,    -1,    84,    85,    86,    87,    88,
      89,    90,    91,    -1,    -1,    -1,    -1,    -1,    -1,    96,
      97,    98,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   131,    -1,    -1,    -1,    -1,    -1,    -1,   116,
     117,    -1,    -1,    -1,    -1,   144,    -1,    -1,    -1,    -1,
      -1,    -1,   129,   130,   131,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   142,   143,   144,   145,   146,
     147,   148,   149,   150,   151,   152,   153,   154,   155,   156,
     157,   158,    -1,    -1,   161,    -1,    -1,    -1,    -1,   168,
      -1,    -1,   169,   170,   171,   172,   173,   174,   175,    -1,
      -1,     3,     4,     5,     6,     7,     8,   184,    10,    11,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    20,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      32,    33,    34,    35,    36,    37,    38,    39,    -1,    41,
      42,    43,    -1,    45,    46,    -1,    -1,    49,    50,    51,
      52,    -1,    54,    -1,    -1,    -1,    -1,    59,    60,    61,
      62,    63,    64,    65,    66,    67,    68,    69,    70,    71,
      72,    73,    74,    75,    76,    77,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    96,    97,    98,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,     3,   116,   117,     6,     7,     8,    -1,
      10,    11,    -1,    -1,    -1,    -1,    -1,   129,   130,   131,
      20,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     142,   143,   144,   145,   146,   147,   148,   149,   150,   151,
     152,   153,   154,   155,   156,   157,   158,    -1,    -1,   161,
      -1,    -1,    -1,    -1,    54,    -1,    -1,   169,   170,   171,
     172,   173,   174,   175,    -1,    -1,     5,    -1,    -1,    -1,
      -1,    -1,   184,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      80,    81,    82,    -1,    84,    85,    86,    87,    88,    89,
      90,    91,    -1,    -1,    -1,    -1,    96,    97,    98,    38,
      39,     3,    -1,    -1,     6,     7,     8,    -1,    10,    11,
      49,    50,    51,    52,    -1,    -1,   116,   117,    20,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   129,
      69,    70,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   143,    -1,   145,    -1,    -1,    -1,   149,
     150,   151,    54,    -1,    -1,    -1,    -1,    -1,   158,    -1,
      -1,   161,    -1,    -1,    -1,    -1,    -1,    -1,   168,   169,
     170,   171,   172,   173,   174,   175,    -1,    -1,    -1,    -1,
      -1,    -1,     3,    -1,   184,     6,     7,     8,    -1,    10,
      11,    -1,   131,    -1,    96,    97,    98,    -1,    -1,    20,
      -1,    -1,    -1,    -1,    -1,   144,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   116,   117,    -1,    -1,     3,    -1,
      -1,     6,     7,     8,    -1,    10,    11,   129,    -1,    -1,
      -1,    -1,    -1,    54,    -1,    20,    -1,    -1,    -1,    -1,
      -1,   143,    -1,   145,    -1,    -1,    -1,   149,   150,   151,
      -1,    -1,    -1,    -1,    -1,    -1,   158,    -1,    -1,   161,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   169,   170,   171,
     172,   173,   174,   175,    -1,    96,    97,    98,    -1,    -1,
       3,    -1,   184,     6,     7,     8,    -1,    10,    11,    -1,
      -1,    -1,    -1,    -1,    -1,   116,   117,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   129,    -1,
      -1,    96,    -1,    98,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   143,    -1,   145,    -1,    -1,    -1,   149,   150,
     151,   116,   117,    -1,    -1,    -1,    -1,   158,    -1,    -1,
     161,    -1,    -1,    -1,   129,    -1,    -1,    -1,    -1,   170,
     171,   172,   173,   174,   175,    -1,    -1,    -1,   143,    -1,
     145,    -1,    -1,   184,   149,   150,   151,    -1,    -1,    -1,
       5,    -1,    -1,    96,    -1,    98,   161,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   170,   171,   172,   173,   174,
     175,    -1,    -1,   116,   117,    -1,    -1,    -1,    -1,   184,
      -1,    -1,    -1,    38,    39,    -1,   129,    -1,    -1,    -1,
      45,    -1,    -1,    -1,    49,    50,    51,    52,    -1,    -1,
     143,    -1,   145,    -1,    -1,    -1,   149,   150,   151,    -1,
      -1,    -1,    -1,    -1,    69,    70,    -1,    -1,   161,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   170,   171,   172,
     173,   174,   175,    -1,    -1,     3,     4,     5,     6,     7,
       8,   184,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    20,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    32,    33,    34,    35,    36,    37,
      38,    39,    -1,    41,    42,    43,   131,    45,    46,    -1,
      -1,    49,    50,    51,    52,    -1,    -1,    -1,    -1,   144,
      -1,    59,    60,    61,    62,    63,    64,    65,    66,    67,
      68,    69,    70,    71,    72,    73,    74,    75,    76,    77,
      -1,    -1,    -1,    -1,    -1,    -1,   171,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    96,    -1,
      98,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   116,   117,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   129,   130,   131,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   142,   143,   144,   145,   146,   147,
     148,   149,   150,   151,   152,   153,   154,   155,   156,   157,
      -1,    -1,    -1,   161,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   170,   171,   172,   173,   174,   175,     3,     4,
       5,     6,     7,     8,    -1,    10,    11,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    32,    33,    34,
      35,    36,    37,    38,    39,    -1,    41,    42,    43,    -1,
      45,    46,    -1,    -1,    49,    50,    51,    52,    -1,    -1,
      -1,    -1,    -1,    -1,    59,    60,    61,    62,    63,    64,
      65,    66,    67,    68,    69,    70,    71,    72,    73,    74,
      75,    76,    77,    -1,    -1,    -1,    -1,    -1,     3,    -1,
      -1,     6,     7,     8,    -1,    10,    11,    -1,    -1,    -1,
      -1,    96,    -1,    98,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   116,   117,    -1,     3,    -1,    -1,     6,     7,     8,
      -1,    10,    11,    -1,   129,   130,   131,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,   142,   143,   144,
     145,   146,    -1,   148,   149,   150,   151,   152,   153,   154,
     155,   156,   157,    -1,    -1,    -1,   161,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   170,   171,   172,   173,   174,
     175,    96,    -1,    98,     3,    -1,    -1,     6,     7,     8,
      -1,    10,    11,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   116,   117,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   129,    -1,    -1,    96,    -1,    98,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   143,    -1,
     145,    -1,    -1,    -1,   149,   150,   151,   116,   117,    -1,
      -1,    -1,    -1,    -1,     3,    -1,   161,     6,     7,     8,
     129,    10,    11,   168,    -1,   170,   171,   172,   173,   174,
     175,    -1,    -1,    -1,   143,    -1,   145,    -1,    -1,    -1,
     149,   150,   151,    -1,    -1,    -1,    -1,    96,    -1,    98,
      -1,     3,   161,    -1,     6,     7,     8,   166,    10,    11,
      -1,   170,   171,   172,   173,   174,   175,   116,   117,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     129,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   143,    -1,   145,    -1,    -1,    -1,
     149,   150,   151,    -1,    -1,    -1,    -1,    96,    -1,    98,
      -1,    -1,   161,   162,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   170,   171,   172,   173,   174,   175,   116,   117,    -1,
      -1,     3,    -1,    -1,     6,     7,     8,    -1,    10,    11,
     129,    -1,    -1,    -1,    96,    -1,    98,    -1,    20,    -1,
      -1,    -1,    -1,    -1,   143,    -1,   145,    -1,    -1,    -1,
     149,   150,   151,    -1,   116,   117,    -1,    -1,    -1,    -1,
      -1,    -1,   161,    -1,    -1,    -1,    -1,   129,    -1,   168,
      -1,   170,   171,   172,   173,   174,   175,    -1,    -1,    -1,
      -1,   143,    -1,   145,    -1,    -1,    -1,   149,   150,   151,
       3,    -1,    -1,     6,     7,     8,    -1,    10,    11,   161,
      -1,    -1,    -1,    -1,   166,    -1,    -1,    -1,   170,   171,
     172,   173,   174,   175,    96,    -1,    98,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   116,   117,    -1,    -1,     3,    -1,
      -1,     6,     7,     8,    -1,    10,    11,   129,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   143,    -1,   145,    -1,    -1,    -1,   149,   150,   151,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   161,
      -1,    -1,    -1,    96,    -1,    98,    -1,    -1,   170,   171,
     172,   173,   174,   175,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,     3,   116,   117,     6,     7,     8,    -1,    10,
      11,    -1,    -1,    -1,    -1,    -1,   129,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     143,    96,   145,    98,    -1,    -1,   149,   150,   151,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   161,    -1,
      -1,   116,   117,    -1,    -1,    -1,    -1,   170,   171,   172,
     173,   174,   175,    -1,   129,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   143,    -1,
     145,    -1,    -1,    -1,   149,   150,   151,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    96,   161,    98,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   170,   171,   172,   173,   174,
     175,    -1,    -1,    -1,    -1,   116,   117,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   129,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   143,    -1,   145,    -1,    -1,    -1,   149,   150,
     151,     3,     4,     5,    -1,    -1,    -1,    -1,    -1,    -1,
     161,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    20,   170,
     171,   172,   173,   174,   175,    -1,    -1,    -1,    -1,    -1,
      32,    33,    34,    35,    36,    37,    38,    39,    -1,    41,
      42,    43,    44,    45,    46,    -1,    -1,    49,    50,    51,
      52,     3,     4,     5,    -1,    -1,    -1,    59,    60,    61,
      62,    63,    64,    65,    66,    67,    68,    69,    70,    71,
      72,    73,    74,    75,    76,    77,    -1,    -1,    -1,    -1,
      32,    33,    34,    35,    36,    37,    38,    39,    -1,    41,
      42,    43,    -1,    45,    46,    -1,    -1,    49,    50,    51,
      52,    -1,    -1,    -1,    -1,    -1,    -1,    59,    60,    61,
      62,    63,    64,    65,    66,    67,    68,    69,    70,    71,
      72,    73,    74,    75,    76,    77,    -1,    -1,   130,   131,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     142,    -1,   144,    -1,   146,   147,   148,    -1,    -1,    -1,
     152,   153,   154,   155,   156,   157,    -1,    -1,    -1,   161,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   171,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   130,   131,
      -1,    -1,    -1,     4,     5,    -1,    -1,    -1,    -1,    -1,
     142,    -1,   144,    -1,   146,    -1,   148,    -1,    -1,    -1,
     152,   153,   154,   155,   156,   157,    -1,    -1,    -1,   161,
     162,    32,    33,    34,    35,    36,    37,    38,    39,   171,
      41,    42,    43,    -1,    45,    46,    -1,    -1,    49,    50,
      51,    52,     4,     5,    -1,    -1,    -1,    -1,    59,    60,
      61,    62,    63,    64,    65,    66,    67,    68,    69,    70,
      71,    72,    73,    74,    75,    76,    77,    -1,    -1,    -1,
      32,    33,    34,    35,    36,    37,    38,    39,    -1,    41,
      42,    43,    -1,    45,    46,    -1,    -1,    49,    50,    51,
      52,    -1,    -1,    -1,    -1,    -1,    -1,    59,    60,    61,
      62,    63,    64,    65,    66,    67,    68,    69,    70,    71,
      72,    73,    74,    75,    76,    77,    -1,    -1,    -1,   130,
     131,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   142,    -1,   144,    -1,   146,    -1,   148,    -1,    -1,
      -1,   152,   153,   154,   155,   156,   157,    -1,    -1,    -1,
     161,   162,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     171,    -1,    -1,    -1,    -1,    -1,    -1,    -1,   130,   131,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
     142,    -1,   144,    -1,   146,    -1,   148,     4,     5,    -1,
     152,   153,   154,   155,   156,   157,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    20,    -1,    -1,   168,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    32,    33,    34,    35,    36,
      37,    38,    39,    -1,    41,    42,    43,    -1,    45,    46,
      -1,    -1,    49,    50,    51,    52,     4,     5,    -1,    -1,
      -1,    -1,    59,    60,    61,    62,    63,    64,    65,    66,
      67,    68,    69,    70,    71,    72,    73,    74,    75,    76,
      77,    -1,    -1,    -1,    32,    33,    34,    35,    36,    37,
      38,    39,    -1,    41,    42,    43,    -1,    45,    46,    -1,
      -1,    49,    50,    51,    52,    53,    -1,    -1,    -1,    -1,
      -1,    59,    60,    61,    62,    63,    64,    65,    66,    67,
      68,    69,    70,    71,    72,    73,    74,    75,    76,    77,
      -1,    -1,    -1,   130,   131,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,   142,    -1,   144,    -1,   146,
     147,   148,    -1,    -1,    -1,   152,   153,   154,   155,   156,
     157,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   130,   131,     4,     5,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,   142,    -1,   144,    -1,   146,    -1,
     148,    -1,    -1,    -1,   152,   153,   154,   155,   156,   157,
      -1,    -1,    32,    33,    34,    35,    36,    37,    38,    39,
      -1,    41,    42,    43,    -1,    45,    46,    -1,    -1,    49,
      50,    51,    52,    -1,     4,     5,    -1,    -1,    -1,    59,
      60,    61,    62,    63,    64,    65,    66,    67,    68,    69,
      70,    71,    72,    73,    74,    75,    76,    77,    -1,    -1,
      -1,    81,    32,    33,    34,    35,    36,    37,    38,    39,
      -1,    41,    42,    43,    -1,    45,    46,    -1,    -1,    49,
      50,    51,    52,    -1,    -1,    -1,    -1,    -1,    -1,    59,
      60,    61,    62,    63,    64,    65,    66,    67,    68,    69,
      70,    71,    72,    73,    74,    75,    76,    77,    -1,    -1,
     130,   131,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   142,    -1,   144,    -1,   146,    -1,   148,    -1,
      -1,    -1,   152,   153,   154,   155,   156,   157,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    40,    -1,    -1,     4,     5,    -1,    -1,    47,    48,
     130,   131,    -1,    -1,    -1,    54,    55,    56,    57,    58,
      -1,    -1,   142,    -1,   144,    -1,   146,    -1,   148,    -1,
      -1,    -1,   152,   153,   154,   155,   156,   157,    38,    39,
      -1,    41,    42,    43,    -1,    45,    46,    -1,    -1,    49,
      50,    51,    52,    92,    93,    94,    95,    -1,    -1,    59,
      60,    61,    62,    63,    64,    65,    66,    67,    68,    69,
      70,    71,    72,    73,    74,    75,    76,    77,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,   132,   133,   134,   135,   136,   137,   138,
     139,   140,   141,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,   131,    -1,    -1,    -1,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   142,    -1,   144,    -1,    -1,    -1,    -1,    -1,
      -1,    -1,   152,   153,   154,   155,   156,   157
};

/* YYSTOS[STATE-NUM] -- The symbol kind of the accessing symbol of
   state STATE-NUM.  */
static const yytype_int16 yystos[] =
{
       0,     3,     4,     5,    20,    32,    33,    34,    35,    36,
      37,    38,    39,    41,    42,    43,    44,    45,    46,    49,
      50,    51,    52,    59,    60,    61,    62,    63,    64,    65,
      66,    67,    68,    69,    70,    71,    72,    73,    74,    75,
      76,    77,   130,   131,   142,   144,   146,   147,   148,   152,
     153,   154,   155,   156,   157,   161,   171,   218,   219,   220,
     223,   224,   225,   226,   227,   229,   232,   235,   239,   240,
     241,   242,   243,   244,   245,   246,   247,   249,   250,   268,
     269,   271,   299,   300,   301,   305,   312,   313,   336,   337,
     272,   273,   274,   275,   337,    56,   337,     6,   161,   168,
     337,   161,   161,   161,   161,   161,   242,   184,   221,   222,
     242,   220,   269,   270,   220,   220,   302,   165,   161,   270,
     270,    40,    47,    48,    54,    55,    56,    57,    58,    92,
      93,    94,    95,   132,   133,   134,   135,   136,   137,   138,
     139,   140,   141,   306,   307,   311,   243,    45,   239,   249,
     251,   184,   271,   300,     0,   270,   163,   166,   161,   276,
      21,    69,   337,     6,     7,     8,    10,    11,    96,    98,
     116,   117,   129,   143,   145,   149,   150,   151,   161,   170,
     171,   172,   173,   174,   175,   186,   187,   188,   191,   195,
     196,   197,   198,   199,   200,   201,   202,   203,   204,   205,
     206,   208,   210,   217,   310,   337,   237,   238,   337,   168,
     217,   219,   320,   217,   195,   210,   212,   214,   320,   214,
     320,   217,   162,   163,   184,   183,   303,   314,   315,   337,
     168,   219,   286,   288,   304,   328,   329,    34,   251,   252,
     162,   248,   319,   337,   161,   161,   161,   161,   161,   161,
     307,   161,    45,   239,   249,   274,   166,     6,    96,   162,
     277,   278,   337,   337,   337,   184,   161,   195,   161,   161,
     195,   195,   195,   195,   195,   161,   161,   214,   320,     9,
      10,    11,   161,   165,   167,   197,   171,   176,   177,   172,
     173,    12,    13,    14,    15,   178,   179,    16,    17,   170,
     180,   181,    18,    19,   182,   162,   163,   236,   270,   237,
     162,   161,   249,   256,   257,   258,   259,   260,   162,   162,
     163,    22,    23,    24,    25,    26,    27,    28,    29,    30,
      31,   183,   213,   162,   163,   162,   162,   162,   162,   222,
     242,   168,   212,   262,   304,   228,   221,    54,    97,   158,
     169,   214,   216,   218,   219,   269,   270,   280,   284,   287,
     289,   290,   321,   322,   323,   327,   330,   335,   337,   288,
     328,   251,   252,    34,   166,   212,   219,   253,   254,   255,
     320,   162,   163,   217,   217,   217,   319,   217,   217,   217,
     161,   162,   278,   184,   320,   320,   320,   320,   212,   162,
     162,   193,   162,   194,   212,   214,   192,   197,   197,   197,
     198,   198,   199,   199,   200,   200,   200,   200,   201,   201,
     202,   203,   204,   207,   209,   211,   238,   169,   183,   316,
     236,   162,   253,   256,   257,   161,   165,   184,    96,   212,
     215,   165,   167,   169,   263,   264,   265,   266,   267,   168,
     184,   184,   161,   184,   214,    80,    81,    82,    84,    85,
      86,    87,    88,    89,    90,    91,   281,   288,   292,   295,
     298,   331,   332,   333,   287,   269,   279,   280,   283,   284,
     218,   323,   280,   284,   324,   325,   326,   164,   212,   212,
     166,   161,   242,   249,   162,   163,   337,   162,   162,   162,
     162,   163,   162,   162,   217,   162,   163,   162,   162,   163,
     168,   197,   337,   162,   163,   166,   337,   205,   206,   214,
     217,   169,   162,   162,   162,   261,   166,   217,   162,   212,
     217,   337,   163,   169,   262,   183,   267,   230,   231,   270,
      96,   334,   184,   217,   285,   161,   161,   337,   184,   184,
     184,   214,   161,   279,   282,   161,   279,   323,   327,   166,
     166,    53,   255,   217,   162,   308,   337,    81,   189,   190,
     320,   169,   263,   194,   164,   253,   166,   184,   166,   169,
     264,   169,   231,   232,   317,   318,   162,   164,   164,   214,
     214,   184,   184,   214,    85,   216,   218,   162,   162,   165,
     167,   164,   162,   163,   164,   163,   169,   210,   162,   262,
     164,   233,   234,   242,   232,   184,   162,   162,   162,   161,
     184,   216,   214,   309,   212,   190,   212,   169,   217,   163,
     184,   164,   293,   294,   296,   214,   216,   184,   166,   337,
     234,   217,   279,   282,   282,   162,   184,   216,    83,   291,
     184,   216,   162,   279,   162,   297,   282,   282
};

/* YYR1[RULE-NUM] -- Symbol kind of the left-hand side of rule RULE-NUM.  */
static const yytype_int16 yyr1[] =
{
       0,   185,   186,   186,   186,   186,   186,   186,   187,   187,
     187,   188,   189,   189,   190,   190,   191,   191,   191,   191,
     192,   191,   193,   191,   191,   191,   191,   191,   191,   194,
     194,   195,   195,   195,   195,   195,   195,   195,   195,   195,
     195,   195,   196,   196,   196,   196,   196,   196,   197,   197,
     198,   198,   198,   198,   199,   199,   199,   200,   200,   200,
     201,   201,   201,   201,   201,   202,   202,   202,   203,   203,
     204,   204,   205,   205,   206,   207,   206,   208,   209,   208,
     210,   211,   210,   212,   212,   213,   213,   213,   213,   213,
     213,   213,   213,   213,   213,   213,   214,   215,   214,   216,
     216,   217,   218,   218,   218,   218,   219,   220,   220,   220,
     220,   220,   220,   221,   221,   222,   222,   223,   224,   224,
     224,   224,   224,   224,   224,   225,   225,   225,   225,   225,
     225,   225,   225,   225,   225,   225,   225,   225,   225,   225,
     225,   225,   225,   225,   225,   225,   225,   225,   225,   225,
     226,   226,   226,   226,   228,   227,   227,   229,   229,   230,
     230,   231,   232,   232,   232,   233,   233,   234,   234,   234,
     234,   235,   235,   235,   236,   236,   237,   237,   238,   239,
     239,   239,   239,   239,   239,   239,   239,   239,   239,   239,
     240,   240,   241,   241,   242,   242,   243,   243,   243,   243,
     244,   244,   244,   245,   245,   245,   245,   246,   246,   247,
     248,   247,   247,   249,   249,   249,   249,   249,   249,   250,
     251,   251,   252,   252,   253,   253,   254,   254,   255,   255,
     256,   256,   256,   257,   257,   257,   258,   258,   259,   259,
     260,   260,   260,   261,   260,   262,   262,   262,   262,   263,
     263,   264,   264,   265,   266,   266,   267,   267,   268,   268,
     269,   269,   270,   270,   271,   272,   272,   273,   273,   274,
     274,   275,   275,   276,   276,   277,   277,   278,   278,   278,
     279,   279,   280,   280,   280,   280,   280,   281,   281,   281,
     282,   283,   284,   284,   285,   284,   286,   287,   288,   288,
     289,   289,   289,   290,   290,   291,   291,   293,   292,   294,
     292,   296,   295,   295,   295,   297,   295,   298,   298,   298,
     298,   298,   299,   299,   300,   300,   300,   302,   301,   303,
     301,   304,   304,   305,   305,   306,   306,   307,   307,   307,
     307,   307,   307,   307,   307,   307,   307,   307,   307,   307,
     307,   307,   307,   307,   307,   307,   307,   307,   307,   308,
     309,   308,   308,   310,   310,   311,   311,   312,   312,   313,
     313,   313,   313,   314,   314,   315,   316,   316,   317,   318,
     318,   319,   319,   320,   320,   321,   322,   323,   323,   325,
     324,   326,   326,   327,   327,   328,   329,   329,   330,   330,
     330,   330,   331,   332,   333,   334,   335,   335,   336,   336,
     337
};

/* YYR2[RULE-NUM] -- Number of symbols on the right-hand side of rule RULE-NUM.  */
static const yytype_int8 yyr2[] =
{
       0,     2,     1,     1,     1,     3,     1,     1,     1,     1,
       1,     6,     1,     3,     3,     3,     1,     4,     3,     4,
       0,     4,     0,     4,     2,     2,     6,     7,     5,     1,
       3,     1,     2,     2,     2,     2,     4,     4,     6,     2,
       2,     2,     1,     1,     1,     1,     1,     1,     1,     4,
       1,     3,     3,     3,     1,     3,     3,     1,     3,     3,
       1,     3,     3,     3,     3,     1,     3,     3,     1,     3,
       1,     3,     1,     3,     1,     0,     4,     1,     0,     4,
       1,     0,     6,     1,     3,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     0,     4,     0,
       1,     1,     2,     3,     1,     1,     1,     1,     2,     1,
       2,     1,     2,     1,     3,     1,     3,     2,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     4,     1,     1,     1,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     2,     4,     1,
       4,     4,     4,     4,     0,     7,     3,     1,     1,     1,
       2,     4,     1,     1,     1,     1,     3,     1,     2,     3,
       0,     5,     6,     2,     0,     1,     1,     3,     3,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     4,     4,     1,     2,     2,     3,     2,     2,
       1,     3,     1,     4,     5,     6,     6,     1,     2,     3,
       0,     5,     4,     1,     5,     2,     6,     2,     3,     1,
       1,     2,     0,     1,     1,     3,     1,     3,     2,     1,
       1,     1,     2,     3,     1,     1,     0,     1,     3,     4,
       2,     3,     3,     0,     5,     1,     2,     3,     4,     2,
       4,     0,     1,     2,     1,     2,     3,     2,     7,     5,
       2,     1,     0,     1,     4,     1,     3,     1,     2,     0,
       1,     1,     3,     2,     3,     1,     2,     1,     1,     1,
       1,     1,     1,     2,     2,     1,     1,     1,     1,     1,
       1,     2,     2,     4,     0,     4,     1,     1,     2,     3,
       1,     1,     2,     2,     3,     2,     0,     0,     7,     0,
       6,     0,     6,     7,     9,     0,     9,     3,     2,     2,
       2,     3,     1,     2,     1,     1,     1,     0,     3,     0,
       4,     1,     2,     0,     1,     1,     2,     4,     1,     1,
       1,     1,     1,     1,     1,     1,     1,     1,     1,     1,
       1,     1,     1,     4,     1,     1,     4,     6,     4,     1,
       0,     4,     4,     1,     1,     1,     4,     1,     1,     1,
       2,     1,     1,     1,     0,     1,     2,     0,     1,     1,
       2,     1,     3,     1,     2,     1,     2,     1,     2,     0,
       2,     2,     1,     1,     2,     3,     1,     2,     1,     1,
       2,     2,     1,     1,     1,     1,     5,     2,     4,     5,
       1
};


enum { YYENOMEM = -2 };

#define yyerrok         (yyerrstatus = 0)
#define yyclearin       (yychar = YYEMPTY)

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYNOMEM         goto yyexhaustedlab


#define YYRECOVERING()  (!!yyerrstatus)

#define YYBACKUP(Token, Value)                                    \
  do                                                              \
    if (yychar == YYEMPTY)                                        \
      {                                                           \
        yychar = (Token);                                         \
        yylval = (Value);                                         \
        YYPOPSTACK (yylen);                                       \
        yystate = *yyssp;                                         \
        goto yybackup;                                            \
      }                                                           \
    else                                                          \
      {                                                           \
        yyerror (YY_("syntax error: cannot back up")); \
        YYERROR;                                                  \
      }                                                           \
  while (0)

/* Backward compatibility with an undocumented macro.
   Use YYerror or YYUNDEF. */
#define YYERRCODE YYUNDEF


/* Enable debugging if requested.  */
#if YYDEBUG

# ifndef YYFPRINTF
#  include <stdio.h> /* INFRINGES ON USER NAME SPACE */
#  define YYFPRINTF fprintf
# endif

# define YYDPRINTF(Args)                        \
do {                                            \
  if (yydebug)                                  \
    YYFPRINTF Args;                             \
} while (0)




# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)                    \
do {                                                                      \
  if (yydebug)                                                            \
    {                                                                     \
      YYFPRINTF (stderr, "%s ", Title);                                   \
      yy_symbol_print (stderr,                                            \
                  Kind, Value); \
      YYFPRINTF (stderr, "\n");                                           \
    }                                                                     \
} while (0)


/*-----------------------------------.
| Print this symbol's value on YYO.  |
`-----------------------------------*/

static void
yy_symbol_value_print (FILE *yyo,
                       yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  FILE *yyoutput = yyo;
  YY_USE (yyoutput);
  if (!yyvaluep)
    return;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/*---------------------------.
| Print this symbol on YYO.  |
`---------------------------*/

static void
yy_symbol_print (FILE *yyo,
                 yysymbol_kind_t yykind, YYSTYPE const * const yyvaluep)
{
  YYFPRINTF (yyo, "%s %s (",
             yykind < YYNTOKENS ? "token" : "nterm", yysymbol_name (yykind));

  yy_symbol_value_print (yyo, yykind, yyvaluep);
  YYFPRINTF (yyo, ")");
}

/*------------------------------------------------------------------.
| yy_stack_print -- Print the state stack from its BOTTOM up to its |
| TOP (included).                                                   |
`------------------------------------------------------------------*/

static void
yy_stack_print (yy_state_t *yybottom, yy_state_t *yytop)
{
  YYFPRINTF (stderr, "Stack now");
  for (; yybottom <= yytop; yybottom++)
    {
      int yybot = *yybottom;
      YYFPRINTF (stderr, " %d", yybot);
    }
  YYFPRINTF (stderr, "\n");
}

# define YY_STACK_PRINT(Bottom, Top)                            \
do {                                                            \
  if (yydebug)                                                  \
    yy_stack_print ((Bottom), (Top));                           \
} while (0)


/*------------------------------------------------.
| Report that the YYRULE is going to be reduced.  |
`------------------------------------------------*/

static void
yy_reduce_print (yy_state_t *yyssp, YYSTYPE *yyvsp,
                 int yyrule)
{
  int yylno = yyrline[yyrule];
  int yynrhs = yyr2[yyrule];
  int yyi;
  YYFPRINTF (stderr, "Reducing stack by rule %d (line %d):\n",
             yyrule - 1, yylno);
  /* The symbols being reduced.  */
  for (yyi = 0; yyi < yynrhs; yyi++)
    {
      YYFPRINTF (stderr, "   $%d = ", yyi + 1);
      yy_symbol_print (stderr,
                       YY_ACCESSING_SYMBOL (+yyssp[yyi + 1 - yynrhs]),
                       &yyvsp[(yyi + 1) - (yynrhs)]);
      YYFPRINTF (stderr, "\n");
    }
}

# define YY_REDUCE_PRINT(Rule)          \
do {                                    \
  if (yydebug)                          \
    yy_reduce_print (yyssp, yyvsp, Rule); \
} while (0)

/* Nonzero means print parse trace.  It is left uninitialized so that
   multiple parsers can coexist.  */
int yydebug;
#else /* !YYDEBUG */
# define YYDPRINTF(Args) ((void) 0)
# define YY_SYMBOL_PRINT(Title, Kind, Value, Location)
# define YY_STACK_PRINT(Bottom, Top)
# define YY_REDUCE_PRINT(Rule)
#endif /* !YYDEBUG */


/* YYINITDEPTH -- initial size of the parser's stacks.  */
#ifndef YYINITDEPTH
# define YYINITDEPTH 200
#endif

/* YYMAXDEPTH -- maximum size the stacks can grow to (effective only
   if the built-in stack extension method is used).

   Do not make this value too large; the results are undefined if
   YYSTACK_ALLOC_MAXIMUM < YYSTACK_BYTES (YYMAXDEPTH)
   evaluated with infinite-precision integer arithmetic.  */

#ifndef YYMAXDEPTH
# define YYMAXDEPTH 10000
#endif






/*-----------------------------------------------.
| Release the memory associated to this symbol.  |
`-----------------------------------------------*/

static void
yydestruct (const char *yymsg,
            yysymbol_kind_t yykind, YYSTYPE *yyvaluep)
{
  YY_USE (yyvaluep);
  if (!yymsg)
    yymsg = "Deleting";
  YY_SYMBOL_PRINT (yymsg, yykind, yyvaluep, yylocationp);

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  YY_USE (yykind);
  YY_IGNORE_MAYBE_UNINITIALIZED_END
}


/* Lookahead token kind.  */
int yychar;

/* The semantic value of the lookahead symbol.  */
YYSTYPE yylval;
/* Number of syntax errors so far.  */
int yynerrs;




/*----------.
| yyparse.  |
`----------*/

int
yyparse (void)
{
    yy_state_fast_t yystate = 0;
    /* Number of tokens to shift before error messages enabled.  */
    int yyerrstatus = 0;

    /* Refer to the stacks through separate pointers, to allow yyoverflow
       to reallocate them elsewhere.  */

    /* Their size.  */
    YYPTRDIFF_T yystacksize = YYINITDEPTH;

    /* The state stack: array, bottom, top.  */
    yy_state_t yyssa[YYINITDEPTH];
    yy_state_t *yyss = yyssa;
    yy_state_t *yyssp = yyss;

    /* The semantic value stack: array, bottom, top.  */
    YYSTYPE yyvsa[YYINITDEPTH];
    YYSTYPE *yyvs = yyvsa;
    YYSTYPE *yyvsp = yyvs;

  int yyn;
  /* The return value of yyparse.  */
  int yyresult;
  /* Lookahead symbol kind.  */
  yysymbol_kind_t yytoken = YYSYMBOL_YYEMPTY;
  /* The variables used to return semantic value and location from the
     action routines.  */
  YYSTYPE yyval;



#define YYPOPSTACK(N)   (yyvsp -= (N), yyssp -= (N))

  /* The number of symbols on the RHS of the reduced rule.
     Keep to zero when no symbol should be popped.  */
  int yylen = 0;

  YYDPRINTF ((stderr, "Starting parse\n"));

  yychar = YYEMPTY; /* Cause a token to be read.  */

  goto yysetstate;


/*------------------------------------------------------------.
| yynewstate -- push a new state, which is found in yystate.  |
`------------------------------------------------------------*/
yynewstate:
  /* In all cases, when you get here, the value and location stacks
     have just been pushed.  So pushing a state here evens the stacks.  */
  yyssp++;


/*--------------------------------------------------------------------.
| yysetstate -- set current state (the top of the stack) to yystate.  |
`--------------------------------------------------------------------*/
yysetstate:
  YYDPRINTF ((stderr, "Entering state %d\n", yystate));
  YY_ASSERT (0 <= yystate && yystate < YYNSTATES);
  YY_IGNORE_USELESS_CAST_BEGIN
  *yyssp = YY_CAST (yy_state_t, yystate);
  YY_IGNORE_USELESS_CAST_END
  YY_STACK_PRINT (yyss, yyssp);

  if (yyss + yystacksize - 1 <= yyssp)
#if !defined yyoverflow && !defined YYSTACK_RELOCATE
    YYNOMEM;
#else
    {
      /* Get the current used size of the three stacks, in elements.  */
      YYPTRDIFF_T yysize = yyssp - yyss + 1;

# if defined yyoverflow
      {
        /* Give user a chance to reallocate the stack.  Use copies of
           these so that the &'s don't force the real ones into
           memory.  */
        yy_state_t *yyss1 = yyss;
        YYSTYPE *yyvs1 = yyvs;

        /* Each stack pointer address is followed by the size of the
           data in use in that stack, in bytes.  This used to be a
           conditional around just the two extra args, but that might
           be undefined if yyoverflow is a macro.  */
        yyoverflow (YY_("memory exhausted"),
                    &yyss1, yysize * YYSIZEOF (*yyssp),
                    &yyvs1, yysize * YYSIZEOF (*yyvsp),
                    &yystacksize);
        yyss = yyss1;
        yyvs = yyvs1;
      }
# else /* defined YYSTACK_RELOCATE */
      /* Extend the stack our own way.  */
      if (YYMAXDEPTH <= yystacksize)
        YYNOMEM;
      yystacksize *= 2;
      if (YYMAXDEPTH < yystacksize)
        yystacksize = YYMAXDEPTH;

      {
        yy_state_t *yyss1 = yyss;
        union yyalloc *yyptr =
          YY_CAST (union yyalloc *,
                   YYSTACK_ALLOC (YY_CAST (YYSIZE_T, YYSTACK_BYTES (yystacksize))));
        if (! yyptr)
          YYNOMEM;
        YYSTACK_RELOCATE (yyss_alloc, yyss);
        YYSTACK_RELOCATE (yyvs_alloc, yyvs);
#  undef YYSTACK_RELOCATE
        if (yyss1 != yyssa)
          YYSTACK_FREE (yyss1);
      }
# endif

      yyssp = yyss + yysize - 1;
      yyvsp = yyvs + yysize - 1;

      YY_IGNORE_USELESS_CAST_BEGIN
      YYDPRINTF ((stderr, "Stack size increased to %ld\n",
                  YY_CAST (long, yystacksize)));
      YY_IGNORE_USELESS_CAST_END

      if (yyss + yystacksize - 1 <= yyssp)
        YYABORT;
    }
#endif /* !defined yyoverflow && !defined YYSTACK_RELOCATE */


  if (yystate == YYFINAL)
    YYACCEPT;

  goto yybackup;


/*-----------.
| yybackup.  |
`-----------*/
yybackup:
  /* Do appropriate processing given the current state.  Read a
     lookahead token if we need one and don't already have one.  */

  /* First try to decide what to do without reference to lookahead token.  */
  yyn = yypact[yystate];
  if (yypact_value_is_default (yyn))
    goto yydefault;

  /* Not known => get a lookahead token if don't already have one.  */

  /* YYCHAR is either empty, or end-of-input, or a valid lookahead.  */
  if (yychar == YYEMPTY)
    {
      YYDPRINTF ((stderr, "Reading a token\n"));
      yychar = yylex ();
    }

  if (yychar <= YYEOF)
    {
      yychar = YYEOF;
      yytoken = YYSYMBOL_YYEOF;
      YYDPRINTF ((stderr, "Now at end of input.\n"));
    }
  else if (yychar == YYerror)
    {
      /* The scanner already issued an error message, process directly
         to error recovery.  But do not keep the error token as
         lookahead, it is too special and may lead us to an endless
         loop in error recovery. */
      yychar = YYUNDEF;
      yytoken = YYSYMBOL_YYerror;
      goto yyerrlab1;
    }
  else
    {
      yytoken = YYTRANSLATE (yychar);
      YY_SYMBOL_PRINT ("Next token is", yytoken, &yylval, &yylloc);
    }

  /* If the proper action on seeing token YYTOKEN is to reduce or to
     detect an error, take that action.  */
  yyn += yytoken;
  if (yyn < 0 || YYLAST < yyn || yycheck[yyn] != yytoken)
    goto yydefault;
  yyn = yytable[yyn];
  if (yyn <= 0)
    {
      if (yytable_value_is_error (yyn))
        goto yyerrlab;
      yyn = -yyn;
      goto yyreduce;
    }

  /* Count tokens shifted since error; after three, turn off error
     status.  */
  if (yyerrstatus)
    yyerrstatus--;

  /* Shift the lookahead token.  */
  YY_SYMBOL_PRINT ("Shifting", yytoken, &yylval, &yylloc);
  yystate = yyn;
  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END

  /* Discard the shifted token.  */
  yychar = YYEMPTY;
  goto yynewstate;


/*-----------------------------------------------------------.
| yydefault -- do the default action for the current state.  |
`-----------------------------------------------------------*/
yydefault:
  yyn = yydefact[yystate];
  if (yyn == 0)
    goto yyerrlab;
  goto yyreduce;


/*-----------------------------.
| yyreduce -- do a reduction.  |
`-----------------------------*/
yyreduce:
  /* yyn is the number of a rule to reduce with.  */
  yylen = yyr2[yyn];

  /* If YYLEN is nonzero, implement the default value of the action:
     '$$ = $1'.

     Otherwise, the following line sets YYVAL to garbage.
     This behavior is undocumented and Bison
     users should not rely upon it.  Assigning to YYVAL
     unconditionally makes the parser a bit smaller, and it avoids a
     GCC warning that YYVAL may be used uninitialized.  */
  yyval = yyvsp[1-yylen];


  YY_REDUCE_PRINT (yyn);
  switch (yyn)
    {
  case 2: /* primary_expression: identifier  */
#line 174 "SDCC.y"
                     { (yyval.asts) = newAst_VALUE (symbolVal ((yyvsp[0].sym))); }
#line 2870 "y.tab.c"
    break;

  case 3: /* primary_expression: CONSTANT  */
#line 175 "SDCC.y"
                     { (yyval.asts) = newAst_VALUE ((yyvsp[0].val)); }
#line 2876 "y.tab.c"
    break;

  case 5: /* primary_expression: '(' expression ')'  */
#line 177 "SDCC.y"
                           { (yyval.asts) = (yyvsp[-1].asts); }
#line 2882 "y.tab.c"
    break;

  case 8: /* predefined_constant: TOKEN_FALSE  */
#line 183 "SDCC.y"
                 { (yyval.asts) = newAst_VALUE (constBoolVal (false, true)); }
#line 2888 "y.tab.c"
    break;

  case 9: /* predefined_constant: TOKEN_TRUE  */
#line 184 "SDCC.y"
                 { (yyval.asts) = newAst_VALUE (constBoolVal (true, true)); }
#line 2894 "y.tab.c"
    break;

  case 10: /* predefined_constant: NULLPTR  */
#line 185 "SDCC.y"
                 { (yyval.asts) = newAst_VALUE (constNullptrVal ()); }
#line 2900 "y.tab.c"
    break;

  case 11: /* generic_selection: GENERIC '(' assignment_expr ',' generic_assoc_list ')'  */
#line 189 "SDCC.y"
                                                            { (yyval.asts) = newNode (GENERIC, (yyvsp[-3].asts), (yyvsp[-1].asts)); }
#line 2906 "y.tab.c"
    break;

  case 12: /* generic_assoc_list: generic_association  */
#line 193 "SDCC.y"
                         { (yyval.asts) = newNode  (GENERIC_ASSOC_LIST, NULL, (yyvsp[0].asts)); }
#line 2912 "y.tab.c"
    break;

  case 13: /* generic_assoc_list: generic_assoc_list ',' generic_association  */
#line 194 "SDCC.y"
                                                { (yyval.asts) = newNode  (GENERIC_ASSOC_LIST, (yyvsp[-2].asts), (yyvsp[0].asts)); }
#line 2918 "y.tab.c"
    break;

  case 14: /* generic_association: type_name ':' assignment_expr  */
#line 198 "SDCC.y"
                                   { (yyval.asts) = newNode  (GENERIC_ASSOCIATION, newAst_LINK((yyvsp[-2].lnk)), (yyvsp[0].asts)); }
#line 2924 "y.tab.c"
    break;

  case 15: /* generic_association: DEFAULT ':' assignment_expr  */
#line 199 "SDCC.y"
                                 { (yyval.asts) = newNode  (GENERIC_ASSOCIATION,NULL,(yyvsp[0].asts)); }
#line 2930 "y.tab.c"
    break;

  case 17: /* postfix_expression: postfix_expression '[' expression ']'  */
#line 204 "SDCC.y"
                                                    { (yyval.asts) = newNode  ('[', (yyvsp[-3].asts), (yyvsp[-1].asts)); }
#line 2936 "y.tab.c"
    break;

  case 18: /* postfix_expression: postfix_expression '(' ')'  */
#line 205 "SDCC.y"
                                              { (yyval.asts) = newNode  (CALL,(yyvsp[-2].asts),NULL);
                                          (yyval.asts)->left->funcName = 1;}
#line 2943 "y.tab.c"
    break;

  case 19: /* postfix_expression: postfix_expression '(' argument_expr_list ')'  */
#line 208 "SDCC.y"
          {
            (yyval.asts) = newNode  (CALL,(yyvsp[-3].asts),(yyvsp[-1].asts)); (yyval.asts)->left->funcName = 1;
          }
#line 2951 "y.tab.c"
    break;

  case 20: /* $@1: %empty  */
#line 211 "SDCC.y"
                            { ignoreTypedefType = 1; }
#line 2957 "y.tab.c"
    break;

  case 21: /* postfix_expression: postfix_expression '.' $@1 identifier  */
#line 212 "SDCC.y"
                      {
                        ignoreTypedefType = 0;
                        (yyvsp[0].sym) = newSymbol((yyvsp[0].sym)->name,NestLevel);
                        (yyvsp[0].sym)->implicit = 1;
                        (yyval.asts) = newNode(PTR_OP,newNode('&',(yyvsp[-3].asts),NULL),newAst_VALUE(symbolVal((yyvsp[0].sym))));
                      }
#line 2968 "y.tab.c"
    break;

  case 22: /* $@2: %empty  */
#line 218 "SDCC.y"
                               { ignoreTypedefType = 1; }
#line 2974 "y.tab.c"
    break;

  case 23: /* postfix_expression: postfix_expression PTR_OP $@2 identifier  */
#line 219 "SDCC.y"
                      {
                        ignoreTypedefType = 0;
                        (yyvsp[0].sym) = newSymbol((yyvsp[0].sym)->name,NestLevel);
                        (yyvsp[0].sym)->implicit = 1;
                        (yyval.asts) = newNode(PTR_OP,(yyvsp[-3].asts),newAst_VALUE(symbolVal((yyvsp[0].sym))));
                      }
#line 2985 "y.tab.c"
    break;

  case 24: /* postfix_expression: postfix_expression INC_OP  */
#line 226 "SDCC.y"
                      { (yyval.asts) = newNode(INC_OP,(yyvsp[-1].asts),NULL);}
#line 2991 "y.tab.c"
    break;

  case 25: /* postfix_expression: postfix_expression DEC_OP  */
#line 228 "SDCC.y"
                      { (yyval.asts) = newNode(DEC_OP,(yyvsp[-1].asts),NULL); }
#line 2997 "y.tab.c"
    break;

  case 26: /* postfix_expression: '(' type_name ')' '{' initializer_list '}'  */
#line 230 "SDCC.y"
                      {
                        /* if (!options.std_c99) */
                          werror(E_COMPOUND_LITERALS_C99);

                        /* TODO: implement compound literals (C99) */
                        (yyval.asts) = newAst_VALUE (valueFromLit (0));
                      }
#line 3009 "y.tab.c"
    break;

  case 27: /* postfix_expression: '(' type_name ')' '{' initializer_list ',' '}'  */
#line 238 "SDCC.y"
     {
       // if (!options.std_c99)
         werror(E_COMPOUND_LITERALS_C99);

       // TODO: implement compound literals (C99)
       (yyval.asts) = newAst_VALUE (valueFromLit (0));
     }
#line 3021 "y.tab.c"
    break;

  case 28: /* postfix_expression: '(' type_name ')' '{' '}'  */
#line 246 "SDCC.y"
     {
       if (!options.std_c2x)
         werror(W_EMPTY_INIT_C2X);
       // if (!options.std_c99)
         werror(E_COMPOUND_LITERALS_C99);

       // TODO: implement compound literals (C99)
       (yyval.asts) = newAst_VALUE (valueFromLit (0));
     }
#line 3035 "y.tab.c"
    break;

  case 30: /* argument_expr_list: assignment_expr ',' argument_expr_list  */
#line 259 "SDCC.y"
                                            { (yyval.asts) = newNode(PARAM,(yyvsp[-2].asts),(yyvsp[0].asts)); }
#line 3041 "y.tab.c"
    break;

  case 32: /* unary_expression: INC_OP unary_expression  */
#line 264 "SDCC.y"
                                    { (yyval.asts) = newNode (INC_OP, NULL, (yyvsp[0].asts)); }
#line 3047 "y.tab.c"
    break;

  case 33: /* unary_expression: DEC_OP unary_expression  */
#line 265 "SDCC.y"
                                    { (yyval.asts) = newNode (DEC_OP, NULL, (yyvsp[0].asts)); }
#line 3053 "y.tab.c"
    break;

  case 34: /* unary_expression: unary_operator cast_expression  */
#line 267 "SDCC.y"
       {
         if ((yyvsp[-1].yyint) == '&' && IS_AST_OP ((yyvsp[0].asts)) && (yyvsp[0].asts)->opval.op == '*' && (yyvsp[0].asts)->right == NULL)
           (yyval.asts) = (yyvsp[0].asts)->left;
         else if ((yyvsp[-1].yyint) == '*' && IS_AST_OP ((yyvsp[0].asts)) && (yyvsp[0].asts)->opval.op == '&' && (yyvsp[0].asts)->right == NULL)
           (yyval.asts) = (yyvsp[0].asts)->left;
         else
           (yyval.asts) = newNode ((yyvsp[-1].yyint), (yyvsp[0].asts), NULL);
       }
#line 3066 "y.tab.c"
    break;

  case 35: /* unary_expression: SIZEOF unary_expression  */
#line 275 "SDCC.y"
                                    { (yyval.asts) = newNode (SIZEOF, NULL, (yyvsp[0].asts)); }
#line 3072 "y.tab.c"
    break;

  case 36: /* unary_expression: SIZEOF '(' type_name ')'  */
#line 276 "SDCC.y"
                                         { (yyval.asts) = newAst_VALUE (sizeofOp ((yyvsp[-1].lnk))); }
#line 3078 "y.tab.c"
    break;

  case 37: /* unary_expression: ALIGNOF '(' type_name ')'  */
#line 277 "SDCC.y"
                                         { (yyval.asts) = newAst_VALUE (alignofOp ((yyvsp[-1].lnk))); }
#line 3084 "y.tab.c"
    break;

  case 38: /* unary_expression: OFFSETOF '(' type_name ',' offsetof_member_designator ')'  */
#line 278 "SDCC.y"
                                                               { (yyval.asts) = offsetofOp((yyvsp[-3].lnk), (yyvsp[-1].asts)); }
#line 3090 "y.tab.c"
    break;

  case 39: /* unary_expression: RLC unary_expression  */
#line 279 "SDCC.y"
                                    { (yyval.asts) = newNode (RLC, (yyvsp[0].asts), NULL); }
#line 3096 "y.tab.c"
    break;

  case 40: /* unary_expression: RRC unary_expression  */
#line 280 "SDCC.y"
                                    { (yyval.asts) = newNode (RRC, (yyvsp[0].asts), NULL); }
#line 3102 "y.tab.c"
    break;

  case 41: /* unary_expression: SWAP unary_expression  */
#line 281 "SDCC.y"
                                    { (yyval.asts) = newNode (SWAP, (yyvsp[0].asts), NULL); }
#line 3108 "y.tab.c"
    break;

  case 42: /* unary_operator: '&'  */
#line 285 "SDCC.y"
            { (yyval.yyint) = '&';}
#line 3114 "y.tab.c"
    break;

  case 43: /* unary_operator: '*'  */
#line 286 "SDCC.y"
            { (yyval.yyint) = '*';}
#line 3120 "y.tab.c"
    break;

  case 44: /* unary_operator: '+'  */
#line 287 "SDCC.y"
            { (yyval.yyint) = '+';}
#line 3126 "y.tab.c"
    break;

  case 45: /* unary_operator: '-'  */
#line 288 "SDCC.y"
            { (yyval.yyint) = '-';}
#line 3132 "y.tab.c"
    break;

  case 46: /* unary_operator: '~'  */
#line 289 "SDCC.y"
            { (yyval.yyint) = '~';}
#line 3138 "y.tab.c"
    break;

  case 47: /* unary_operator: '!'  */
#line 290 "SDCC.y"
            { (yyval.yyint) = '!';}
#line 3144 "y.tab.c"
    break;

  case 49: /* cast_expression: '(' type_name ')' cast_expression  */
#line 295 "SDCC.y"
                                       { (yyval.asts) = newNode(CAST,newAst_LINK((yyvsp[-2].lnk)),(yyvsp[0].asts)); }
#line 3150 "y.tab.c"
    break;

  case 51: /* multiplicative_expression: multiplicative_expression '*' cast_expression  */
#line 300 "SDCC.y"
                                                   { (yyval.asts) = newNode('*',(yyvsp[-2].asts),(yyvsp[0].asts));}
#line 3156 "y.tab.c"
    break;

  case 52: /* multiplicative_expression: multiplicative_expression '/' cast_expression  */
#line 301 "SDCC.y"
                                                   { (yyval.asts) = newNode('/',(yyvsp[-2].asts),(yyvsp[0].asts));}
#line 3162 "y.tab.c"
    break;

  case 53: /* multiplicative_expression: multiplicative_expression '%' cast_expression  */
#line 302 "SDCC.y"
                                                   { (yyval.asts) = newNode('%',(yyvsp[-2].asts),(yyvsp[0].asts));}
#line 3168 "y.tab.c"
    break;

  case 55: /* additive_expression: additive_expression '+' multiplicative_expression  */
#line 307 "SDCC.y"
                                                       { (yyval.asts)=newNode('+',(yyvsp[-2].asts),(yyvsp[0].asts));}
#line 3174 "y.tab.c"
    break;

  case 56: /* additive_expression: additive_expression '-' multiplicative_expression  */
#line 308 "SDCC.y"
                                                       { (yyval.asts)=newNode('-',(yyvsp[-2].asts),(yyvsp[0].asts));}
#line 3180 "y.tab.c"
    break;

  case 58: /* shift_expression: shift_expression LEFT_OP additive_expression  */
#line 313 "SDCC.y"
                                                   { (yyval.asts) = newNode(LEFT_OP,(yyvsp[-2].asts),(yyvsp[0].asts)); }
#line 3186 "y.tab.c"
    break;

  case 59: /* shift_expression: shift_expression RIGHT_OP additive_expression  */
#line 314 "SDCC.y"
                                                   { (yyval.asts) = newNode(RIGHT_OP,(yyvsp[-2].asts),(yyvsp[0].asts)); }
#line 3192 "y.tab.c"
    break;

  case 61: /* relational_expression: relational_expression '<' shift_expression  */
#line 319 "SDCC.y"
                                                  { (yyval.asts) = newNode('<',  (yyvsp[-2].asts),(yyvsp[0].asts));}
#line 3198 "y.tab.c"
    break;

  case 62: /* relational_expression: relational_expression '>' shift_expression  */
#line 320 "SDCC.y"
                                                  { (yyval.asts) = newNode('>',  (yyvsp[-2].asts),(yyvsp[0].asts));}
#line 3204 "y.tab.c"
    break;

  case 63: /* relational_expression: relational_expression LE_OP shift_expression  */
#line 321 "SDCC.y"
                                                  { (yyval.asts) = newNode(LE_OP,(yyvsp[-2].asts),(yyvsp[0].asts));}
#line 3210 "y.tab.c"
    break;

  case 64: /* relational_expression: relational_expression GE_OP shift_expression  */
#line 322 "SDCC.y"
                                                  { (yyval.asts) = newNode(GE_OP,(yyvsp[-2].asts),(yyvsp[0].asts));}
#line 3216 "y.tab.c"
    break;

  case 66: /* equality_expression: equality_expression EQ_OP relational_expression  */
#line 327 "SDCC.y"
                                                     { (yyval.asts) = newNode(EQ_OP,(yyvsp[-2].asts),(yyvsp[0].asts));}
#line 3222 "y.tab.c"
    break;

  case 67: /* equality_expression: equality_expression NE_OP relational_expression  */
#line 328 "SDCC.y"
                                                     { (yyval.asts) = newNode(NE_OP,(yyvsp[-2].asts),(yyvsp[0].asts));}
#line 3228 "y.tab.c"
    break;

  case 69: /* and_expression: and_expression '&' equality_expression  */
#line 333 "SDCC.y"
                                             { (yyval.asts) = newNode('&',(yyvsp[-2].asts),(yyvsp[0].asts));}
#line 3234 "y.tab.c"
    break;

  case 71: /* exclusive_or_expression: exclusive_or_expression '^' and_expression  */
#line 338 "SDCC.y"
                                                { (yyval.asts) = newNode('^',(yyvsp[-2].asts),(yyvsp[0].asts));}
#line 3240 "y.tab.c"
    break;

  case 73: /* inclusive_or_expression: inclusive_or_expression '|' exclusive_or_expression  */
#line 343 "SDCC.y"
                                                         { (yyval.asts) = newNode('|',(yyvsp[-2].asts),(yyvsp[0].asts));}
#line 3246 "y.tab.c"
    break;

  case 75: /* $@3: %empty  */
#line 348 "SDCC.y"
                             { seqPointNo++;}
#line 3252 "y.tab.c"
    break;

  case 76: /* logical_and_expr: logical_and_expr AND_OP $@3 inclusive_or_expression  */
#line 349 "SDCC.y"
                                 { (yyval.asts) = newNode(AND_OP,(yyvsp[-3].asts),(yyvsp[0].asts));}
#line 3258 "y.tab.c"
    break;

  case 78: /* $@4: %empty  */
#line 354 "SDCC.y"
                           { seqPointNo++;}
#line 3264 "y.tab.c"
    break;

  case 79: /* logical_or_expr: logical_or_expr OR_OP $@4 logical_and_expr  */
#line 355 "SDCC.y"
                                 { (yyval.asts) = newNode(OR_OP,(yyvsp[-3].asts),(yyvsp[0].asts)); }
#line 3270 "y.tab.c"
    break;

  case 81: /* $@5: %empty  */
#line 360 "SDCC.y"
                         { seqPointNo++;}
#line 3276 "y.tab.c"
    break;

  case 82: /* conditional_expr: logical_or_expr '?' $@5 expression ':' conditional_expr  */
#line 361 "SDCC.y"
                     {
                        (yyval.asts) = newNode(':',(yyvsp[-2].asts),(yyvsp[0].asts));
                        (yyval.asts) = newNode('?',(yyvsp[-5].asts),(yyval.asts));
                     }
#line 3285 "y.tab.c"
    break;

  case 84: /* assignment_expr: unary_expression assignment_operator assignment_expr  */
#line 370 "SDCC.y"
                     {

                             switch ((yyvsp[-1].yyint)) {
                             case '=':
                                     (yyval.asts) = newNode((yyvsp[-1].yyint),(yyvsp[-2].asts),(yyvsp[0].asts));
                                     break;
                             case MUL_ASSIGN:
                                     (yyval.asts) = createRMW((yyvsp[-2].asts), '*', (yyvsp[0].asts));
                                     break;
                             case DIV_ASSIGN:
                                     (yyval.asts) = createRMW((yyvsp[-2].asts), '/', (yyvsp[0].asts));
                                     break;
                             case MOD_ASSIGN:
                                     (yyval.asts) = createRMW((yyvsp[-2].asts), '%', (yyvsp[0].asts));
                                     break;
                             case ADD_ASSIGN:
                                     (yyval.asts) = createRMW((yyvsp[-2].asts), '+', (yyvsp[0].asts));
                                     break;
                             case SUB_ASSIGN:
                                     (yyval.asts) = createRMW((yyvsp[-2].asts), '-', (yyvsp[0].asts));
                                     break;
                             case LEFT_ASSIGN:
                                     (yyval.asts) = createRMW((yyvsp[-2].asts), LEFT_OP, (yyvsp[0].asts));
                                     break;
                             case RIGHT_ASSIGN:
                                     (yyval.asts) = createRMW((yyvsp[-2].asts), RIGHT_OP, (yyvsp[0].asts));
                                     break;
                             case AND_ASSIGN:
                                     (yyval.asts) = createRMW((yyvsp[-2].asts), '&', (yyvsp[0].asts));
                                     break;
                             case XOR_ASSIGN:
                                     (yyval.asts) = createRMW((yyvsp[-2].asts), '^', (yyvsp[0].asts));
                                     break;
                             case OR_ASSIGN:
                                     (yyval.asts) = createRMW((yyvsp[-2].asts), '|', (yyvsp[0].asts));
                                     break;
                             default :
                                     (yyval.asts) = NULL;
                             }

                     }
#line 3331 "y.tab.c"
    break;

  case 85: /* assignment_operator: '='  */
#line 414 "SDCC.y"
                     { (yyval.yyint) = '=';}
#line 3337 "y.tab.c"
    break;

  case 97: /* $@6: %empty  */
#line 429 "SDCC.y"
                    { seqPointNo++;}
#line 3343 "y.tab.c"
    break;

  case 98: /* expression: expression ',' $@6 assignment_expr  */
#line 429 "SDCC.y"
                                                     { (yyval.asts) = newNode(',',(yyvsp[-3].asts),(yyvsp[0].asts));}
#line 3349 "y.tab.c"
    break;

  case 99: /* expression_opt: %empty  */
#line 433 "SDCC.y"
                     { (yyval.asts) = NULL; seqPointNo++; }
#line 3355 "y.tab.c"
    break;

  case 100: /* expression_opt: expression  */
#line 434 "SDCC.y"
                     { (yyval.asts) = (yyvsp[0].asts); seqPointNo++; }
#line 3361 "y.tab.c"
    break;

  case 102: /* declaration: declaration_specifiers ';'  */
#line 445 "SDCC.y"
      {
         /* Special case: if incomplete struct/union declared without name, */
         /* make sure an incomplete type for it exists in the current scope */
         if (IS_STRUCT((yyvsp[-1].lnk)))
           {
             structdef *sdef = SPEC_STRUCT((yyvsp[-1].lnk));
             structdef *osdef;
             osdef = findSymWithBlock (StructTab, sdef->tagsym, currBlockno, NestLevel);
             if (osdef && osdef->block != currBlockno)
               {
                 sdef = newStruct(osdef->tagsym->name);
                 sdef->level = NestLevel;
                 sdef->block = currBlockno;
                 sdef->tagsym = newSymbol (osdef->tagsym->name, NestLevel);
                 addSym (StructTab, sdef, sdef->tag, sdef->level, currBlockno, 0);
                 uselessDecl = false;
               }
           }
         if (uselessDecl)
           werror(W_USELESS_DECL);
         uselessDecl = true;
         (yyval.sym) = NULL;
      }
#line 3389 "y.tab.c"
    break;

  case 103: /* declaration: declaration_specifiers init_declarator_list ';'  */
#line 469 "SDCC.y"
      {
         /* add the specifier list to the id */
         symbol *sym , *sym1;

         bool autocandidate = options.std_c2x && IS_SPEC ((yyvsp[-2].lnk)) && SPEC_SCLS((yyvsp[-2].lnk)) == S_AUTO;

         for (sym1 = sym = reverseSyms((yyvsp[-1].sym));sym != NULL;sym = sym->next) {
             sym_link *lnk = copyLinkChain((yyvsp[-2].lnk));
             sym_link *l0 = NULL, *l1 = NULL, *l2 = NULL;
             /* check illegal declaration */
             for (l0 = sym->type; l0 != NULL; l0 = l0->next)
               if (IS_PTR (l0))
                 break;
             /* check if creating instances of structs with flexible arrays */
             for (l1 = lnk; l1 != NULL; l1 = l1->next)
               if (IS_STRUCT (l1) && SPEC_STRUCT (l1)->b_flexArrayMember)
                 break;
             if (!options.std_c99 && l0 == NULL && l1 != NULL && SPEC_EXTR((yyvsp[-2].lnk)) != 1)
               werror (W_FLEXARRAY_INSTRUCT, sym->name);
             /* check if creating instances of function type */
             for (l1 = lnk; l1 != NULL; l1 = l1->next)
               if (IS_FUNC (l1))
                 break;
             for (l2 = lnk; l2 != NULL; l2 = l2->next)
               if (IS_PTR (l2))
                 break;
             if (l0 == NULL && l2 == NULL && l1 != NULL)
               werrorfl(sym->fileDef, sym->lineDef, E_TYPE_IS_FUNCTION, sym->name);
             if (autocandidate && !sym->type && sym->ival && sym->ival->type == INIT_NODE) // C2X auto type inference
               {
                 sym->type = sym->etype = typeofOp (sym->ival->init.node);
                 SPEC_SCLS (lnk) = 0;
               }
             /* do the pointer stuff */
             pointerTypes(sym->type,lnk);
             addDecl (sym,0,lnk);
         }

         uselessDecl = true;
         (yyval.sym) = sym1;
      }
#line 3435 "y.tab.c"
    break;

  case 104: /* declaration: static_assert_declaration  */
#line 511 "SDCC.y"
      {
         (yyval.sym) = NULL;
      }
#line 3443 "y.tab.c"
    break;

  case 105: /* declaration: attribute_declaration  */
#line 515 "SDCC.y"
      {
         (yyval.sym) = NULL;
      }
#line 3451 "y.tab.c"
    break;

  case 106: /* declaration_specifiers: declaration_specifiers_  */
#line 520 "SDCC.y"
                                                 { (yyval.lnk) = finalizeSpec((yyvsp[0].lnk)); }
#line 3457 "y.tab.c"
    break;

  case 107: /* declaration_specifiers_: storage_class_specifier  */
#line 523 "SDCC.y"
                                                     { (yyval.lnk) = (yyvsp[0].lnk); }
#line 3463 "y.tab.c"
    break;

  case 108: /* declaration_specifiers_: storage_class_specifier declaration_specifiers_  */
#line 524 "SDCC.y"
                                                     {
     /* if the decl $2 is not a specifier */
     /* find the spec and replace it      */
     (yyval.lnk) = mergeDeclSpec((yyvsp[-1].lnk), (yyvsp[0].lnk), "storage_class_specifier declaration_specifiers - skipped");
   }
#line 3473 "y.tab.c"
    break;

  case 109: /* declaration_specifiers_: type_specifier_qualifier  */
#line 529 "SDCC.y"
                                                    { (yyval.lnk) = (yyvsp[0].lnk); }
#line 3479 "y.tab.c"
    break;

  case 110: /* declaration_specifiers_: type_specifier_qualifier declaration_specifiers_  */
#line 530 "SDCC.y"
                                                      {
     /* if the decl $2 is not a specifier */
     /* find the spec and replace it      */
     (yyval.lnk) = mergeDeclSpec((yyvsp[-1].lnk), (yyvsp[0].lnk), "type_specifier_qualifier declaration_specifiers - skipped");
   }
#line 3489 "y.tab.c"
    break;

  case 111: /* declaration_specifiers_: function_specifier  */
#line 535 "SDCC.y"
                                                    { (yyval.lnk) = (yyvsp[0].lnk); }
#line 3495 "y.tab.c"
    break;

  case 112: /* declaration_specifiers_: function_specifier declaration_specifiers_  */
#line 536 "SDCC.y"
                                                    {
     /* if the decl $2 is not a specifier */
     /* find the spec and replace it      */
     (yyval.lnk) = mergeDeclSpec((yyvsp[-1].lnk), (yyvsp[0].lnk), "function_specifier declaration_specifiers - skipped");
   }
#line 3505 "y.tab.c"
    break;

  case 114: /* init_declarator_list: init_declarator_list ',' init_declarator  */
#line 545 "SDCC.y"
                                                   { (yyvsp[0].sym)->next = (yyvsp[-2].sym); (yyval.sym) = (yyvsp[0].sym);}
#line 3511 "y.tab.c"
    break;

  case 115: /* init_declarator: declarator  */
#line 549 "SDCC.y"
                                 { (yyvsp[0].sym)->ival = NULL; }
#line 3517 "y.tab.c"
    break;

  case 116: /* init_declarator: declarator '=' initializer  */
#line 550 "SDCC.y"
                                 { (yyvsp[-2].sym)->ival = (yyvsp[0].ilist); seqPointNo++; }
#line 3523 "y.tab.c"
    break;

  case 118: /* storage_class_specifier: TYPEDEF  */
#line 558 "SDCC.y"
               {
                  (yyval.lnk) = newLink (SPECIFIER);
                  SPEC_TYPEDEF((yyval.lnk)) = 1;
               }
#line 3532 "y.tab.c"
    break;

  case 119: /* storage_class_specifier: EXTERN  */
#line 562 "SDCC.y"
               {
                  (yyval.lnk) = newLink(SPECIFIER);
                  SPEC_EXTR((yyval.lnk)) = 1;
               }
#line 3541 "y.tab.c"
    break;

  case 120: /* storage_class_specifier: STATIC  */
#line 566 "SDCC.y"
               {
                  (yyval.lnk) = newLink (SPECIFIER);
                  SPEC_STAT((yyval.lnk)) = 1;
               }
#line 3550 "y.tab.c"
    break;

  case 121: /* storage_class_specifier: THREAD_LOCAL  */
#line 571 "SDCC.y"
               {
                  (yyval.lnk) = 0;
                  werror(E_THREAD_LOCAL);
               }
#line 3559 "y.tab.c"
    break;

  case 122: /* storage_class_specifier: AUTO  */
#line 575 "SDCC.y"
               {
                  (yyval.lnk) = newLink (SPECIFIER);
                  SPEC_SCLS((yyval.lnk)) = S_AUTO;
               }
#line 3568 "y.tab.c"
    break;

  case 123: /* storage_class_specifier: REGISTER  */
#line 579 "SDCC.y"
               {
                  (yyval.lnk) = newLink (SPECIFIER);
                  SPEC_SCLS((yyval.lnk)) = S_REGISTER;
               }
#line 3577 "y.tab.c"
    break;

  case 124: /* storage_class_specifier: CONSTEXPR  */
#line 583 "SDCC.y"
               {
                  (yyval.lnk) = newLink (SPECIFIER);
                  werror (E_CONSTEXPR);
               }
#line 3586 "y.tab.c"
    break;

  case 125: /* type_specifier: SD_VOID  */
#line 590 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  SPEC_NOUN((yyval.lnk)) = V_VOID;
                  ignoreTypedefType = 1;
               }
#line 3596 "y.tab.c"
    break;

  case 126: /* type_specifier: SD_CHAR  */
#line 595 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  SPEC_NOUN((yyval.lnk)) = V_CHAR;
                  ignoreTypedefType = 1;
               }
#line 3606 "y.tab.c"
    break;

  case 127: /* type_specifier: SD_SHORT  */
#line 600 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  SPEC_SHORT((yyval.lnk)) = 1;
                  ignoreTypedefType = 1;
               }
#line 3616 "y.tab.c"
    break;

  case 128: /* type_specifier: SD_INT  */
#line 605 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  SPEC_NOUN((yyval.lnk)) = V_INT;
                  ignoreTypedefType = 1;
               }
#line 3626 "y.tab.c"
    break;

  case 129: /* type_specifier: SD_LONG  */
#line 610 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  SPEC_LONG((yyval.lnk)) = 1;
                  ignoreTypedefType = 1;
               }
#line 3636 "y.tab.c"
    break;

  case 130: /* type_specifier: SD_FLOAT  */
#line 615 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  SPEC_NOUN((yyval.lnk)) = V_FLOAT;
                  ignoreTypedefType = 1;
               }
#line 3646 "y.tab.c"
    break;

  case 131: /* type_specifier: DOUBLE  */
#line 620 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  SPEC_NOUN((yyval.lnk)) = V_DOUBLE;
                  ignoreTypedefType = 1;
               }
#line 3656 "y.tab.c"
    break;

  case 132: /* type_specifier: SIGNED  */
#line 625 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  (yyval.lnk)->select.s.b_signed = 1;
                  ignoreTypedefType = 1;
               }
#line 3666 "y.tab.c"
    break;

  case 133: /* type_specifier: UNSIGNED  */
#line 630 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  SPEC_USIGN((yyval.lnk)) = 1;
                  ignoreTypedefType = 1;
               }
#line 3676 "y.tab.c"
    break;

  case 134: /* type_specifier: SD_BITINT '(' constant_expr ')'  */
#line 635 "SDCC.y"
                                      {
                  (yyval.lnk)=newLink(SPECIFIER);
                  SPEC_NOUN((yyval.lnk)) = V_BITINT;
                  SPEC_BITINTWIDTH((yyval.lnk)) = (unsigned int) ulFromVal(constExprValue((yyvsp[-1].asts), true));
                  ignoreTypedefType = 1;
               }
#line 3687 "y.tab.c"
    break;

  case 135: /* type_specifier: SD_BOOL  */
#line 641 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  SPEC_NOUN((yyval.lnk)) = V_BOOL;
                  ignoreTypedefType = 1;
               }
#line 3697 "y.tab.c"
    break;

  case 136: /* type_specifier: COMPLEX  */
#line 646 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  werror (E_COMPLEX_UNSUPPORTED);
               }
#line 3706 "y.tab.c"
    break;

  case 137: /* type_specifier: IMAGINARY  */
#line 650 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  werror (E_COMPLEX_UNSUPPORTED);
               }
#line 3715 "y.tab.c"
    break;

  case 138: /* type_specifier: DECIMAL32  */
#line 654 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  werror (E_DECIMAL_FLOAT_UNSUPPORTED);
               }
#line 3724 "y.tab.c"
    break;

  case 139: /* type_specifier: DECIMAL64  */
#line 658 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  werror (E_DECIMAL_FLOAT_UNSUPPORTED);
               }
#line 3733 "y.tab.c"
    break;

  case 140: /* type_specifier: DECIMAL128  */
#line 662 "SDCC.y"
                {
                  (yyval.lnk)=newLink(SPECIFIER);
                  werror (E_DECIMAL_FLOAT_UNSUPPORTED);
               }
#line 3742 "y.tab.c"
    break;

  case 141: /* type_specifier: struct_or_union_specifier  */
#line 666 "SDCC.y"
                                {
                                   uselessDecl = false;
                                   (yyval.lnk) = (yyvsp[0].lnk);
                                   ignoreTypedefType = 1;
                                }
#line 3752 "y.tab.c"
    break;

  case 142: /* type_specifier: enum_specifier  */
#line 671 "SDCC.y"
                        {
                           cenum = NULL;
                           uselessDecl = false;
                           ignoreTypedefType = 1;
                           (yyval.lnk) = (yyvsp[0].lnk);
                        }
#line 3763 "y.tab.c"
    break;

  case 143: /* type_specifier: TYPE_NAME  */
#line 678 "SDCC.y"
         {
            symbol *sym;
            sym_link *p;
            sym = findSym(TypedefTab,NULL,(yyvsp[0].yychar));
            (yyval.lnk) = p = copyLinkChain(sym ? sym->type : NULL);
            SPEC_TYPEDEF(getSpec(p)) = 0;
            ignoreTypedefType = 1;
         }
#line 3776 "y.tab.c"
    break;

  case 144: /* type_specifier: typeof_specifier  */
#line 687 "SDCC.y"
     {
       (yyval.lnk) = (yyvsp[0].lnk);
     }
#line 3784 "y.tab.c"
    break;

  case 145: /* type_specifier: FIXED16X16  */
#line 690 "SDCC.y"
                {
                  (yyval.lnk)=newLink(SPECIFIER);
                  SPEC_NOUN((yyval.lnk)) = V_FIXED16X16;
                  ignoreTypedefType = 1;
               }
#line 3794 "y.tab.c"
    break;

  case 146: /* type_specifier: BIT  */
#line 695 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  SPEC_NOUN((yyval.lnk)) = V_BIT;
                  SPEC_SCLS((yyval.lnk)) = S_BIT;
                  SPEC_BLEN((yyval.lnk)) = 1;
                  SPEC_BSTR((yyval.lnk)) = 0;
                  ignoreTypedefType = 1;
               }
#line 3807 "y.tab.c"
    break;

  case 147: /* type_specifier: AT CONSTANT  */
#line 703 "SDCC.y"
                 {
                  (yyval.lnk)=newLink(SPECIFIER);
                  /* add this to the storage class specifier  */
                  SPEC_ABSA((yyval.lnk)) = 1;   /* set the absolute addr flag */
                  /* now get the abs addr from value */
                  SPEC_ADDR((yyval.lnk)) = (unsigned int) ulFromVal(constExprValue(newAst_VALUE ((yyvsp[0].val)), true));
               }
#line 3819 "y.tab.c"
    break;

  case 148: /* type_specifier: AT '(' constant_expr ')'  */
#line 710 "SDCC.y"
                              {
                  (yyval.lnk)=newLink(SPECIFIER);
                  /* add this to the storage class specifier  */
                  SPEC_ABSA((yyval.lnk)) = 1;   /* set the absolute addr flag */
                  /* now get the abs addr from value */
                  SPEC_ADDR((yyval.lnk)) = (unsigned int) ulFromVal(constExprValue((yyvsp[-1].asts), true));
               }
#line 3831 "y.tab.c"
    break;

  case 150: /* typeof_specifier: TYPEOF '(' expression ')'  */
#line 723 "SDCC.y"
     {
       (yyval.lnk) = typeofOp ((yyvsp[-1].asts));
       wassert ((yyval.lnk));
     }
#line 3840 "y.tab.c"
    break;

  case 151: /* typeof_specifier: TYPEOF '(' type_name ')'  */
#line 728 "SDCC.y"
     {
       checkTypeSanity ((yyvsp[-1].lnk), "(typeof)");
       (yyval.lnk) = (yyvsp[-1].lnk);
     }
#line 3849 "y.tab.c"
    break;

  case 152: /* typeof_specifier: TYPEOF_UNQUAL '(' expression ')'  */
#line 733 "SDCC.y"
     {
       (yyval.lnk) = typeofOp ((yyvsp[-1].asts));
       wassert ((yyval.lnk));
       wassert (IS_SPEC ((yyval.lnk)));
       SPEC_CONST ((yyval.lnk)) = 0;
       SPEC_RESTRICT ((yyval.lnk)) = 0;
       SPEC_VOLATILE ((yyval.lnk)) = 0;
       SPEC_ATOMIC ((yyval.lnk)) = 0;
       SPEC_ADDRSPACE ((yyval.lnk)) = 0;
     }
#line 3864 "y.tab.c"
    break;

  case 153: /* typeof_specifier: TYPEOF_UNQUAL '(' type_name ')'  */
#line 744 "SDCC.y"
     {
       checkTypeSanity ((yyvsp[-1].lnk), "(typeof_unqual)");
       (yyval.lnk) = (yyvsp[-1].lnk);
       wassert (IS_SPEC ((yyval.lnk)));
       SPEC_CONST ((yyval.lnk)) = 0;
       SPEC_RESTRICT ((yyval.lnk)) = 0;
       SPEC_VOLATILE ((yyval.lnk)) = 0;
       SPEC_ATOMIC ((yyval.lnk)) = 0;
       SPEC_ADDRSPACE ((yyval.lnk)) = 0;
     }
#line 3879 "y.tab.c"
    break;

  case 154: /* $@7: %empty  */
#line 757 "SDCC.y"
        {
          structdef *sdef;

          if (!(yyvsp[0].sdef)->tagsym)
            {
              /* no tag given, so new struct def for current scope */
              addSym (StructTab, (yyvsp[0].sdef), (yyvsp[0].sdef)->tag, (yyvsp[0].sdef)->level, currBlockno, 0);
            }
          else
            {
              sdef = findSymWithBlock (StructTab, (yyvsp[0].sdef)->tagsym, currBlockno, NestLevel);
              if (sdef)
                {
                  /* Error if a complete type already defined in this scope */
                  if (sdef->block == currBlockno)
                    {
                      if (sdef->fields)
                        (yyvsp[0].sdef)->redefinition = true;
                      else // We are completing an incomplete type
                        (yyvsp[0].sdef) = sdef;
                    }
                  else
                    {
                      /* There is an existing struct def in an outer scope. */
                      /* Create new struct def for current scope */
                      addSym (StructTab, (yyvsp[0].sdef), (yyvsp[0].sdef)->tag, (yyvsp[0].sdef)->level, currBlockno, 0);
                    }
                }
              else
               {
                 /* There is no existing struct def at all. */
                 /* Create new struct def for current scope */
                 addSym (StructTab, (yyvsp[0].sdef), (yyvsp[0].sdef)->tag, (yyvsp[0].sdef)->level, currBlockno, 0);
               }
            }

          if (!(yyvsp[0].sdef)->type)
            {
              (yyvsp[0].sdef)->type = (yyvsp[-2].yyint);
            }
          else
            {
              if ((yyvsp[0].sdef)->type != (yyvsp[-2].yyint))
                  werror(E_BAD_TAG, (yyvsp[0].sdef)->tag, (yyvsp[-2].yyint)==STRUCT ? "struct" : "union");
            }
        }
#line 3930 "y.tab.c"
    break;

  case 155: /* struct_or_union_specifier: struct_or_union attribute_specifier_sequence_opt opt_stag $@7 '{' member_declaration_list '}'  */
#line 804 "SDCC.y"
        {
          structdef *sdef;
          symbol *sym, *dsym;

          // check for errors in structure members
          for (sym=(yyvsp[-1].sym); sym; sym=sym->next)
            {
              if (IS_ABSOLUTE(sym->etype))
                {
                  werrorfl(sym->fileDef, sym->lineDef, E_NOT_ALLOWED, "'at'");
                  SPEC_ABSA(sym->etype) = 0;
                }
              if (IS_SPEC(sym->etype) && SPEC_SCLS(sym->etype))
                {
                  werrorfl(sym->fileDef, sym->lineDef, E_NOT_ALLOWED, "storage class");
                  printTypeChainRaw (sym->type, NULL);
                  SPEC_SCLS(sym->etype) = 0;
                }
              for (dsym=sym->next; dsym; dsym=dsym->next)
                {
                  if (*dsym->name && strcmp(sym->name, dsym->name)==0)
                    {
                      werrorfl(sym->fileDef, sym->lineDef, E_DUPLICATE_MEMBER,
                               (yyvsp[-6].yyint)==STRUCT ? "struct" : "union", sym->name);
                      werrorfl(dsym->fileDef, dsym->lineDef, E_PREVIOUS_DEF);
                    }
                }
            }

          /* Create a structdef   */
          (yyvsp[-4].sdef)->fields = reverseSyms((yyvsp[-1].sym));        /* link the fields */
          (yyvsp[-4].sdef)->size = compStructSize((yyvsp[-6].yyint), (yyvsp[-4].sdef));   /* update size of  */
          promoteAnonStructs ((yyvsp[-6].yyint), (yyvsp[-4].sdef));

          if ((yyvsp[-4].sdef)->redefinition) // Since C2X, multiple definitions for struct / union are allowed, if they are compatible and have the same tags. The current standard draft N3047 allows redeclarations of unions to have a different order of the members. We don't. The rule in N3047 is now considered a mistake by many, and will hopefully be changed to the SDCC behaviour via a national body comment for the final version of the standard.
            {
              sdef = findSymWithBlock (StructTab, (yyvsp[-4].sdef)->tagsym, currBlockno, NestLevel);
              bool compatible = options.std_c2x && sdef->tagsym && (yyvsp[-4].sdef)->tagsym && !strcmp (sdef->tagsym->name, (yyvsp[-4].sdef)->tagsym->name);
              for (symbol *fieldsym1 = sdef->fields, *fieldsym2 = (yyvsp[-4].sdef)->fields; compatible; fieldsym1 = fieldsym1->next, fieldsym2 = fieldsym2->next)
                {
                  if (!fieldsym1 && !fieldsym2)
                    break;
                  if (!fieldsym1 || !fieldsym2)
                    compatible = false;
                  else if (strcmp (fieldsym1->name, fieldsym2->name))
                    compatible = false;
                  else if (compareType (fieldsym1->type, fieldsym2->type, true) <= 0)
                    compatible = false;
               }
              if (!compatible)
                {
                  werror(E_STRUCT_REDEF_INCOMPATIBLE, (yyvsp[-4].sdef)->tag);
                  werrorfl(sdef->tagsym->fileDef, sdef->tagsym->lineDef, E_PREVIOUS_DEF);
                }
            }
          else
            sdef = (yyvsp[-4].sdef);

          /* Create the specifier */
          (yyval.lnk) = newLink (SPECIFIER);
          SPEC_NOUN((yyval.lnk)) = V_STRUCT;
          SPEC_STRUCT((yyval.lnk))= sdef;
        }
#line 3998 "y.tab.c"
    break;

  case 156: /* struct_or_union_specifier: struct_or_union attribute_specifier_sequence_opt stag  */
#line 868 "SDCC.y"
        {
          structdef *sdef;

          sdef = findSymWithBlock (StructTab, (yyvsp[0].sdef)->tagsym, currBlockno, NestLevel);

          if (sdef)
            (yyvsp[0].sdef) = sdef;
          else
            {
              /* new struct def for current scope */
              addSym (StructTab, (yyvsp[0].sdef), (yyvsp[0].sdef)->tag, (yyvsp[0].sdef)->level, currBlockno, 0);
            }
          (yyval.lnk) = newLink(SPECIFIER);
          SPEC_NOUN((yyval.lnk)) = V_STRUCT;
          SPEC_STRUCT((yyval.lnk)) = (yyvsp[0].sdef);

          if (!(yyvsp[0].sdef)->type)
            {
              (yyvsp[0].sdef)->type = (yyvsp[-2].yyint);
            }
          else
            {
              if ((yyvsp[0].sdef)->type != (yyvsp[-2].yyint))
                  werror(E_BAD_TAG, (yyvsp[0].sdef)->tag, (yyvsp[-2].yyint)==STRUCT ? "struct" : "union");
            }
        }
#line 4029 "y.tab.c"
    break;

  case 157: /* struct_or_union: STRUCT  */
#line 897 "SDCC.y"
                     { (yyval.yyint) = STRUCT; ignoreTypedefType = 1; }
#line 4035 "y.tab.c"
    break;

  case 158: /* struct_or_union: UNION  */
#line 898 "SDCC.y"
                     { (yyval.yyint) = UNION; ignoreTypedefType = 1; }
#line 4041 "y.tab.c"
    break;

  case 160: /* member_declaration_list: member_declaration_list member_declaration  */
#line 904 "SDCC.y"
        {
          symbol *sym = (yyvsp[0].sym);

          /* go to the end of the chain */
          while (sym->next) sym = sym->next;
          sym->next = (yyvsp[-1].sym);

           (yyval.sym) = (yyvsp[0].sym);
        }
#line 4055 "y.tab.c"
    break;

  case 161: /* member_declaration: attribute_specifier_sequence_opt specifier_qualifier_list member_declarator_list ';'  */
#line 917 "SDCC.y"
        {
          /* add this type to all the symbols */
          symbol *sym;
          for ( sym = (yyvsp[-1].sym); sym != NULL; sym = sym->next )
            {
              sym_link *btype = copyLinkChain((yyvsp[-2].lnk));

              pointerTypes(sym->type, btype);
              if (!sym->type)
                {
                  sym->type = btype;
                  sym->etype = getSpec(sym->type);
                }
              else
                  addDecl (sym, 0, btype);
              /* make sure the type is complete and sane */
              checkTypeSanity(sym->etype, sym->name);
            }
          ignoreTypedefType = 0;
          (yyval.sym) = (yyvsp[-1].sym);
        }
#line 4081 "y.tab.c"
    break;

  case 162: /* type_specifier_qualifier: type_specifier  */
#line 941 "SDCC.y"
                         { (yyval.lnk) = (yyvsp[0].lnk); }
#line 4087 "y.tab.c"
    break;

  case 163: /* type_specifier_qualifier: type_qualifier  */
#line 942 "SDCC.y"
                         { (yyval.lnk) = (yyvsp[0].lnk); }
#line 4093 "y.tab.c"
    break;

  case 164: /* type_specifier_qualifier: alignment_specifier  */
#line 943 "SDCC.y"
                         { (yyval.lnk) = (yyvsp[0].lnk); }
#line 4099 "y.tab.c"
    break;

  case 166: /* member_declarator_list: member_declarator_list ',' member_declarator  */
#line 949 "SDCC.y"
        {
          (yyvsp[0].sym)->next  = (yyvsp[-2].sym);
          (yyval.sym) = (yyvsp[0].sym);
        }
#line 4108 "y.tab.c"
    break;

  case 168: /* member_declarator: ':' constant_expr  */
#line 958 "SDCC.y"
        {
          unsigned int bitsize;
          (yyval.sym) = newSymbol (genSymName(NestLevel), NestLevel);
          bitsize = (unsigned int) ulFromVal(constExprValue((yyvsp[0].asts), true));
          if (!bitsize)
              bitsize = BITVAR_PAD;
          (yyval.sym)->bitVar = bitsize;
          (yyval.sym)->bitUnnamed = 1;
        }
#line 4122 "y.tab.c"
    break;

  case 169: /* member_declarator: declarator ':' constant_expr  */
#line 968 "SDCC.y"
        {
          unsigned int bitsize;
          bitsize = (unsigned int) ulFromVal(constExprValue((yyvsp[0].asts), true));

          if (!bitsize)
            {
              (yyval.sym) = newSymbol (genSymName(NestLevel), NestLevel);
              (yyval.sym)->bitVar = BITVAR_PAD;
              werror(W_BITFLD_NAMED);
            }
          else
              (yyvsp[-2].sym)->bitVar = bitsize;
        }
#line 4140 "y.tab.c"
    break;

  case 170: /* member_declarator: %empty  */
#line 981 "SDCC.y"
     { (yyval.sym) = newSymbol ("", NestLevel); }
#line 4146 "y.tab.c"
    break;

  case 171: /* enum_specifier: ENUM '{' enumerator_list enum_comma_opt '}'  */
#line 986 "SDCC.y"
        {
          (yyval.lnk) = newEnumType ((yyvsp[-2].sym));
          SPEC_SCLS(getSpec((yyval.lnk))) = 0;
        }
#line 4155 "y.tab.c"
    break;

  case 172: /* enum_specifier: ENUM identifier '{' enumerator_list enum_comma_opt '}'  */
#line 991 "SDCC.y"
        {
          symbol *csym;
          sym_link *enumtype;

          enumtype = newEnumType ((yyvsp[-2].sym));
          SPEC_SCLS(getSpec(enumtype)) = 0;
          (yyvsp[-4].sym)->type = enumtype;

          csym = findSymWithLevel(enumTab, (yyvsp[-4].sym));
          if ((csym && csym->level == (yyvsp[-4].sym)->level))
            {
              if (!options.std_c2x || compareType (csym->type, (yyvsp[-4].sym)->type, true) <= 0)
                {
                  werrorfl((yyvsp[-4].sym)->fileDef, (yyvsp[-4].sym)->lineDef, E_DUPLICATE_TYPEDEF, csym->name);
                  werrorfl(csym->fileDef, csym->lineDef, E_PREVIOUS_DEF);
                }
            }

          /* add this to the enumerator table */
          if (!csym)
              addSym (enumTab, (yyvsp[-4].sym), (yyvsp[-4].sym)->name, (yyvsp[-4].sym)->level, (yyvsp[-4].sym)->block, 0);
          (yyval.lnk) = copyLinkChain(enumtype);
        }
#line 4183 "y.tab.c"
    break;

  case 173: /* enum_specifier: ENUM identifier  */
#line 1015 "SDCC.y"
        {
          symbol *csym;

          /* check the enumerator table */
          if ((csym = findSymWithLevel(enumTab, (yyvsp[0].sym))))
              (yyval.lnk) = copyLinkChain(csym->type);
          else
            {
              (yyval.lnk) = newLink(SPECIFIER);
              SPEC_NOUN((yyval.lnk)) = V_INT;
            }
        }
#line 4200 "y.tab.c"
    break;

  case 174: /* enum_comma_opt: %empty  */
#line 1031 "SDCC.y"
     {
       (yyval.lnk) = NULL;
     }
#line 4208 "y.tab.c"
    break;

  case 175: /* enum_comma_opt: ','  */
#line 1035 "SDCC.y"
     {
       if (!options.std_c99)
         werror (E_ENUM_COMMA_C99);
       (yyval.lnk) = NULL;
     }
#line 4218 "y.tab.c"
    break;

  case 177: /* enumerator_list: enumerator_list ',' enumerator  */
#line 1044 "SDCC.y"
        {
          (yyvsp[0].sym)->next = (yyvsp[-2].sym);
          (yyval.sym) = (yyvsp[0].sym);
        }
#line 4227 "y.tab.c"
    break;

  case 178: /* enumerator: identifier attribute_specifier_sequence_opt opt_assign_expr  */
#line 1052 "SDCC.y"
        {
          symbol *sym;

          (yyvsp[-2].sym)->type = copyLinkChain ((yyvsp[0].val)->type);
          (yyvsp[-2].sym)->etype = getSpec ((yyvsp[-2].sym)->type);
          SPEC_ENUM ((yyvsp[-2].sym)->etype) = 1;
          (yyval.sym) = (yyvsp[-2].sym);

          // check if the symbol at the same level already exists
          if ((sym = findSymWithLevel (SymbolTab, (yyvsp[-2].sym))) && sym->level == (yyvsp[-2].sym)->level)
            {
              // C2X allows redefinitions of enumeration constants with the same value as part of a redeclaration of the same enumerated type.
              if (!options.std_c2x || ullFromVal (valFromType (sym->type)) != ullFromVal (valFromType ((yyvsp[-2].sym)->type)))
                {
                  werrorfl ((yyvsp[-2].sym)->fileDef, (yyvsp[-2].sym)->lineDef, E_DUPLICATE_MEMBER, "enum", (yyvsp[-2].sym)->name);
                  werrorfl (sym->fileDef, sym->lineDef, E_PREVIOUS_DEF);
                }
            }
          else
            {
              // do this now, so we can use it for the next enums in the list
              addSymChain (&(yyvsp[-2].sym));
            }
        }
#line 4256 "y.tab.c"
    break;

  case 179: /* type_qualifier: SD_CONST  */
#line 1079 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  SPEC_CONST((yyval.lnk)) = 1;
               }
#line 4265 "y.tab.c"
    break;

  case 180: /* type_qualifier: RESTRICT  */
#line 1083 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  SPEC_RESTRICT((yyval.lnk)) = 1;
               }
#line 4274 "y.tab.c"
    break;

  case 181: /* type_qualifier: VOLATILE  */
#line 1087 "SDCC.y"
               {
                  (yyval.lnk)=newLink(SPECIFIER);
                  SPEC_VOLATILE((yyval.lnk)) = 1;
               }
#line 4283 "y.tab.c"
    break;

  case 182: /* type_qualifier: ATOMIC  */
#line 1091 "SDCC.y"
             {
                  (yyval.lnk)=newLink(SPECIFIER);
                  werror (E_ATOMIC_UNSUPPORTED);
               }
#line 4292 "y.tab.c"
    break;

  case 183: /* type_qualifier: ADDRSPACE_NAME  */
#line 1095 "SDCC.y"
                    {
                  (yyval.lnk)=newLink(SPECIFIER);
                  SPEC_ADDRSPACE((yyval.lnk)) = findSym (AddrspaceTab, 0, (yyvsp[0].yychar));
               }
#line 4301 "y.tab.c"
    break;

  case 184: /* type_qualifier: XDATA  */
#line 1099 "SDCC.y"
               {
                  (yyval.lnk) = newLink (SPECIFIER);
                  SPEC_SCLS((yyval.lnk)) = S_XDATA;
               }
#line 4310 "y.tab.c"
    break;

  case 185: /* type_qualifier: CODE  */
#line 1103 "SDCC.y"
               {
                  (yyval.lnk) = newLink (SPECIFIER);
                  SPEC_SCLS((yyval.lnk)) = S_CODE;
               }
#line 4319 "y.tab.c"
    break;

  case 186: /* type_qualifier: EEPROM  */
#line 1107 "SDCC.y"
               {
                  (yyval.lnk) = newLink (SPECIFIER);
                  SPEC_SCLS((yyval.lnk)) = S_EEPROM;
               }
#line 4328 "y.tab.c"
    break;

  case 187: /* type_qualifier: DATA  */
#line 1111 "SDCC.y"
               {
                  (yyval.lnk) = newLink (SPECIFIER);
                  SPEC_SCLS((yyval.lnk)) = S_DATA;
               }
#line 4337 "y.tab.c"
    break;

  case 188: /* type_qualifier: IDATA  */
#line 1115 "SDCC.y"
               {
                  (yyval.lnk) = newLink (SPECIFIER);
                  SPEC_SCLS((yyval.lnk)) = S_IDATA;
               }
#line 4346 "y.tab.c"
    break;

  case 189: /* type_qualifier: PDATA  */
#line 1119 "SDCC.y"
               {
                  (yyval.lnk) = newLink (SPECIFIER);
                  SPEC_SCLS((yyval.lnk)) = S_PDATA;
               }
#line 4355 "y.tab.c"
    break;

  case 190: /* function_specifier: INLINE  */
#line 1126 "SDCC.y"
               {
                  (yyval.lnk) = newLink (SPECIFIER);
                  SPEC_INLINE((yyval.lnk)) = 1;
               }
#line 4364 "y.tab.c"
    break;

  case 191: /* function_specifier: NORETURN  */
#line 1130 "SDCC.y"
               {
                  (yyval.lnk) = newLink (SPECIFIER);
                  SPEC_NORETURN((yyval.lnk)) = 1;
               }
#line 4373 "y.tab.c"
    break;

  case 192: /* alignment_specifier: ALIGNAS '(' type_name ')'  */
#line 1138 "SDCC.y"
              {
                 checkTypeSanity ((yyvsp[-1].lnk), "(_Alignas)");
                 (yyval.lnk) = newLink (SPECIFIER);
                 SPEC_ALIGNAS((yyval.lnk)) = 1;
              }
#line 4383 "y.tab.c"
    break;

  case 193: /* alignment_specifier: ALIGNAS '(' constant_expr ')'  */
#line 1144 "SDCC.y"
              {
                 value *val = constExprValue ((yyvsp[-1].asts), true);
                 (yyval.lnk) = newLink (SPECIFIER);
                 SPEC_ALIGNAS((yyval.lnk)) = 0;
                 if (!val)
                   werror (E_CONST_EXPECTED);
                 else if (ulFromVal (val) == 0 || isPowerOf2 (ulFromVal (val)) && ulFromVal (val) <= port->mem.maxextalign)
                   SPEC_ALIGNAS((yyval.lnk)) = ulFromVal(val);
                 else
                   werror (E_ALIGNAS, ulFromVal(val));
              }
#line 4399 "y.tab.c"
    break;

  case 194: /* declarator: direct_declarator  */
#line 1158 "SDCC.y"
                                              { (yyval.sym) = (yyvsp[0].sym); }
#line 4405 "y.tab.c"
    break;

  case 195: /* declarator: pointer direct_declarator  */
#line 1160 "SDCC.y"
         {
             addDecl ((yyvsp[0].sym),0,reverseLink((yyvsp[-1].lnk)));
             (yyval.sym) = (yyvsp[0].sym);
         }
#line 4414 "y.tab.c"
    break;

  case 197: /* direct_declarator: '(' declarator ')'  */
#line 1168 "SDCC.y"
                            { (yyval.sym) = (yyvsp[-1].sym); }
#line 4420 "y.tab.c"
    break;

  case 201: /* declarator2: '(' declarator ')'  */
#line 1175 "SDCC.y"
                            { (yyval.sym) = (yyvsp[-1].sym); }
#line 4426 "y.tab.c"
    break;

  case 203: /* array_declarator: direct_declarator '[' type_qualifier_list_opt ']'  */
#line 1181 "SDCC.y"
     {
       sym_link *p, *n;

       p = newLink (DECLARATOR);
       DCL_TYPE(p) = ARRAY;
       DCL_ELEM(p) = 0;

       if ((yyvsp[-1].lnk))
         {
           if (!options.std_c99)
             werror (E_QUALIFIED_ARRAY_PARAM_C99);

           DCL_PTR_CONST(p) = SPEC_CONST ((yyvsp[-1].lnk));
           DCL_PTR_RESTRICT(p) = SPEC_RESTRICT ((yyvsp[-1].lnk));
           DCL_PTR_VOLATILE(p) = SPEC_VOLATILE ((yyvsp[-1].lnk));
           DCL_PTR_ADDRSPACE(p) = SPEC_ADDRSPACE ((yyvsp[-1].lnk));
           addDecl((yyvsp[-3].sym),0,p);
           n = newLink (SPECIFIER);
           SPEC_NEEDSPAR(n) = 1;
           addDecl((yyvsp[-3].sym),0,n);
         }
       else
         addDecl((yyvsp[-3].sym),0,p);
     }
#line 4455 "y.tab.c"
    break;

  case 204: /* array_declarator: direct_declarator '[' type_qualifier_list_opt assignment_expr ']'  */
#line 1206 "SDCC.y"
     {
       sym_link *p, *n;

       p = newLink (DECLARATOR);
       DCL_TYPE(p) = ARRAY;
       DCL_ELEM_AST (p) = (yyvsp[-1].asts);

       if ((yyvsp[-2].lnk))
         {
           if (!options.std_c99)
             werror (E_QUALIFIED_ARRAY_PARAM_C99);
           DCL_PTR_CONST(p) = SPEC_CONST ((yyvsp[-2].lnk));
           DCL_PTR_RESTRICT(p) = SPEC_RESTRICT ((yyvsp[-2].lnk));
           DCL_PTR_VOLATILE(p) = SPEC_VOLATILE ((yyvsp[-2].lnk));
           DCL_PTR_ADDRSPACE(p) = SPEC_ADDRSPACE ((yyvsp[-2].lnk));
           addDecl((yyvsp[-4].sym), 0, p);
           n = newLink (SPECIFIER);
           SPEC_NEEDSPAR(n) = 1;
           addDecl((yyvsp[-4].sym),0,n);
         }
       else
         addDecl((yyvsp[-4].sym), 0, p);
     }
#line 4483 "y.tab.c"
    break;

  case 205: /* array_declarator: direct_declarator '[' STATIC type_qualifier_list_opt assignment_expr ']'  */
#line 1230 "SDCC.y"
     {
       sym_link *p, *n;

       if (!options.std_c99)
         werror (E_STATIC_ARRAY_PARAM_C99);

       p = newLink (DECLARATOR);
       DCL_TYPE(p) = ARRAY;
       DCL_ELEM_AST (p) = (yyvsp[-1].asts);

       if ((yyvsp[-2].lnk))
         {
           if (!options.std_c99)
             werror (E_QUALIFIED_ARRAY_PARAM_C99);
           DCL_PTR_CONST(p) = SPEC_CONST ((yyvsp[-2].lnk));
           DCL_PTR_RESTRICT(p) = SPEC_RESTRICT ((yyvsp[-2].lnk));
           DCL_PTR_VOLATILE(p) = SPEC_VOLATILE ((yyvsp[-2].lnk));
           DCL_PTR_ADDRSPACE(p) = SPEC_ADDRSPACE ((yyvsp[-2].lnk));
         }
       addDecl((yyvsp[-5].sym), 0, p);
       n = newLink (SPECIFIER);
       SPEC_NEEDSPAR(n) = 1;
       addDecl((yyvsp[-5].sym),0,n);
     }
#line 4512 "y.tab.c"
    break;

  case 206: /* array_declarator: direct_declarator '[' type_qualifier_list STATIC assignment_expr ']'  */
#line 1255 "SDCC.y"
     {
       sym_link *p, *n;

       if (!options.std_c99)
         {
           werror (E_QUALIFIED_ARRAY_PARAM_C99);
           werror (E_STATIC_ARRAY_PARAM_C99);
         }

       p = newLink (DECLARATOR);
       DCL_TYPE(p) = ARRAY;
       DCL_ELEM_AST (p) = (yyvsp[-1].asts);

       DCL_PTR_CONST(p) = SPEC_CONST ((yyvsp[-3].lnk));
       DCL_PTR_RESTRICT(p) = SPEC_RESTRICT ((yyvsp[-3].lnk));
       DCL_PTR_VOLATILE(p) = SPEC_VOLATILE ((yyvsp[-3].lnk));
       DCL_PTR_ADDRSPACE(p) = SPEC_ADDRSPACE ((yyvsp[-3].lnk));
       addDecl((yyvsp[-5].sym), 0, p);
       n = newLink (SPECIFIER);
       SPEC_NEEDSPAR(n) = 1;
       addDecl((yyvsp[-5].sym),0,n);
     }
#line 4539 "y.tab.c"
    break;

  case 207: /* declarator2_function_attributes: function_declarator  */
#line 1280 "SDCC.y"
                                         { (yyval.sym) = (yyvsp[0].sym); }
#line 4545 "y.tab.c"
    break;

  case 208: /* declarator2_function_attributes: function_declarator function_attributes  */
#line 1281 "SDCC.y"
                                              {
           // copy the functionAttributes (not the args and hasVargs !!)
           sym_link *funcType=(yyvsp[-1].sym)->type;

           while (funcType && !IS_FUNC(funcType))
             funcType = funcType->next;

           if (!funcType)
             werror (E_FUNC_ATTR);
           else
             {
               struct value *args = FUNC_ARGS(funcType);
               unsigned hasVargs = FUNC_HASVARARGS(funcType);
               bool noprototype = FUNC_NOPROTOTYPE(funcType);

               memcpy (&funcType->funcAttrs, &(yyvsp[0].lnk)->funcAttrs,
                   sizeof((yyvsp[0].lnk)->funcAttrs));

               FUNC_ARGS(funcType)=args;
               FUNC_HASVARARGS(funcType)=hasVargs;
               FUNC_NOPROTOTYPE(funcType)=noprototype;

               // just to be sure
               memset (&(yyvsp[0].lnk)->funcAttrs, 0,
                   sizeof((yyvsp[0].lnk)->funcAttrs));

               addDecl ((yyvsp[-1].sym),0,(yyvsp[0].lnk));
             }
   }
#line 4579 "y.tab.c"
    break;

  case 209: /* function_declarator: declarator2 '(' ')'  */
#line 1314 "SDCC.y"
     {
       addDecl ((yyvsp[-2].sym), FUNCTION, NULL);

       // Up to C17, this was a function declarator without a prototype.
       if (!options.std_c2x)
         {
           FUNC_NOPROTOTYPE((yyvsp[-2].sym)->type) = true;
           if (!options.lessPedantic)
             werror (W_FUNCDECL_WITH_NO_PROTOTYPE);
         }
     }
#line 4595 "y.tab.c"
    break;

  case 210: /* $@8: %empty  */
#line 1326 "SDCC.y"
        {
          NestLevel += LEVEL_UNIT;
          STACK_PUSH(blockNum, currBlockno);
          btree_add_child(currBlockno, ++blockNo);
          currBlockno = blockNo;
          seqPointNo++; /* not a true sequence point, but helps resolve scope */
        }
#line 4607 "y.tab.c"
    break;

  case 211: /* function_declarator: declarator2 '(' $@8 parameter_type_list ')'  */
#line 1334 "SDCC.y"
        {
          sym_link *funcType;

          bool is_fptr = IS_FUNC((yyvsp[-4].sym)->type); // Already a function, must be a function pointer.

          addDecl ((yyvsp[-4].sym), FUNCTION, NULL);
          funcType = (yyvsp[-4].sym)->type;

          // For a function pointer, the parameter list here is for the returned type.
          if (is_fptr)
            funcType = funcType->next;

          while (funcType && !IS_FUNC(funcType))
            funcType = funcType->next;

          wassert (funcType);

          FUNC_HASVARARGS(funcType) = IS_VARG((yyvsp[-1].val));
          FUNC_ARGS(funcType) = reverseVal((yyvsp[-1].val));

          FUNC_SDCCCALL(funcType) = -1;

          /* nest level was incremented to take care of the parms  */
          leaveBlockScope (currBlockno);
          NestLevel -= LEVEL_UNIT;
          currBlockno = STACK_POP(blockNum);
          seqPointNo++; /* not a true sequence point, but helps resolve scope */

          // if this was a pointer (to a function)
          if (!IS_FUNC((yyvsp[-4].sym)->type))
              cleanUpLevel(SymbolTab, NestLevel + LEVEL_UNIT);

          (yyval.sym) = (yyvsp[-4].sym);
        }
#line 4646 "y.tab.c"
    break;

  case 212: /* function_declarator: declarator2 '(' identifier_list ')'  */
#line 1369 "SDCC.y"
        {
          werror(E_OLD_STYLE,(yyvsp[-3].sym)->name);
          
          addDecl ((yyvsp[-3].sym), FUNCTION, NULL);
          
          (yyval.sym) = (yyvsp[-3].sym);
        }
#line 4658 "y.tab.c"
    break;

  case 213: /* pointer: unqualified_pointer  */
#line 1379 "SDCC.y"
                         { (yyval.lnk) = (yyvsp[0].lnk);}
#line 4664 "y.tab.c"
    break;

  case 214: /* pointer: unqualified_pointer AT '(' constant_expr ')'  */
#line 1381 "SDCC.y"
         {
             sym_link *n = newLink(SPECIFIER);
             /* add this to the storage class specifier  */
             SPEC_ABSA(n) = 1;   /* set the absolute addr flag */
             /* now get the abs addr from value */
             SPEC_ADDR(n) = (unsigned int) ulFromVal(constExprValue((yyvsp[-1].asts),true));
             n->next = (yyvsp[-4].lnk);
             (yyval.lnk) = n;
         }
#line 4678 "y.tab.c"
    break;

  case 215: /* pointer: unqualified_pointer type_qualifier_list  */
#line 1391 "SDCC.y"
         {
             (yyval.lnk) = (yyvsp[-1].lnk);
             if (IS_SPEC((yyvsp[0].lnk))) {
                 DCL_TSPEC((yyvsp[-1].lnk)) = (yyvsp[0].lnk);
                 DCL_PTR_CONST((yyvsp[-1].lnk)) = SPEC_CONST((yyvsp[0].lnk));
                 DCL_PTR_VOLATILE((yyvsp[-1].lnk)) = SPEC_VOLATILE((yyvsp[0].lnk));
                 DCL_PTR_RESTRICT((yyvsp[-1].lnk)) = SPEC_RESTRICT((yyvsp[0].lnk));
                 DCL_PTR_ADDRSPACE((yyvsp[-1].lnk)) = SPEC_ADDRSPACE((yyvsp[0].lnk));
             }
             else
                 werror (W_PTR_TYPE_INVALID);
         }
#line 4695 "y.tab.c"
    break;

  case 216: /* pointer: unqualified_pointer type_qualifier_list AT '(' constant_expr ')'  */
#line 1404 "SDCC.y"
         {
             if (IS_SPEC((yyvsp[-4].lnk))) {
                 DCL_TSPEC((yyvsp[-5].lnk)) = (yyvsp[-4].lnk);
                 DCL_PTR_CONST((yyvsp[-5].lnk)) = SPEC_CONST((yyvsp[-4].lnk));
                 DCL_PTR_VOLATILE((yyvsp[-5].lnk)) = SPEC_VOLATILE((yyvsp[-4].lnk));
                 DCL_PTR_RESTRICT((yyvsp[-5].lnk)) = SPEC_RESTRICT((yyvsp[-4].lnk));
                 DCL_PTR_ADDRSPACE((yyvsp[-5].lnk)) = SPEC_ADDRSPACE((yyvsp[-4].lnk));
             }
             else
                 werror (W_PTR_TYPE_INVALID);

             sym_link *n = newLink(SPECIFIER);
             /* add this to the storage class specifier  */
             SPEC_ABSA(n) = 1;   /* set the absolute addr flag */
             /* now get the abs addr from value */
             SPEC_ADDR(n) = (unsigned int) ulFromVal(constExprValue((yyvsp[-1].asts),true));
             n->next = (yyvsp[-5].lnk);
             (yyval.lnk) = n;
         }
#line 4719 "y.tab.c"
    break;

  case 217: /* pointer: unqualified_pointer pointer  */
#line 1424 "SDCC.y"
         {
             (yyval.lnk) = (yyvsp[-1].lnk);
             (yyval.lnk)->next = (yyvsp[0].lnk);
             DCL_TYPE((yyvsp[0].lnk))=port->unqualified_pointer;
         }
#line 4729 "y.tab.c"
    break;

  case 218: /* pointer: unqualified_pointer type_qualifier_list pointer  */
#line 1430 "SDCC.y"
         {
             (yyval.lnk) = (yyvsp[-2].lnk);
             if (IS_SPEC((yyvsp[-1].lnk)) && DCL_TYPE((yyvsp[0].lnk)) == UPOINTER) {
                 DCL_PTR_CONST((yyvsp[-2].lnk)) = SPEC_CONST((yyvsp[-1].lnk));
                 DCL_PTR_VOLATILE((yyvsp[-2].lnk)) = SPEC_VOLATILE((yyvsp[-1].lnk));
                 DCL_PTR_RESTRICT((yyvsp[-2].lnk)) = SPEC_RESTRICT((yyvsp[-1].lnk));
                 DCL_PTR_ADDRSPACE((yyvsp[-2].lnk)) = SPEC_ADDRSPACE((yyvsp[-1].lnk));
                 switch (SPEC_SCLS((yyvsp[-1].lnk))) {
                 case S_XDATA:
                     DCL_TYPE((yyvsp[0].lnk)) = FPOINTER;
                     break;
                 case S_IDATA:
                     DCL_TYPE((yyvsp[0].lnk)) = IPOINTER;
                     break;
                 case S_PDATA:
                     DCL_TYPE((yyvsp[0].lnk)) = PPOINTER;
                     break;
                 case S_DATA:
                     DCL_TYPE((yyvsp[0].lnk)) = POINTER;
                     break;
                 case S_CODE:
                     DCL_TYPE((yyvsp[0].lnk)) = CPOINTER;
                     break;
                 case S_EEPROM:
                     DCL_TYPE((yyvsp[0].lnk)) = EEPPOINTER;
                     break;
                 default:
                   // this could be just "constant"
                   // werror(W_PTR_TYPE_INVALID);
                     ;
                 }
             }
             else
                 werror (W_PTR_TYPE_INVALID);
             (yyval.lnk)->next = (yyvsp[0].lnk);
         }
#line 4770 "y.tab.c"
    break;

  case 219: /* unqualified_pointer: '*'  */
#line 1470 "SDCC.y"
      {
        (yyval.lnk) = newLink(DECLARATOR);
        DCL_TYPE((yyval.lnk))=UPOINTER;
      }
#line 4779 "y.tab.c"
    break;

  case 221: /* type_qualifier_list: type_qualifier_list type_qualifier  */
#line 1479 "SDCC.y"
               {
                 (yyval.lnk) = mergeDeclSpec((yyvsp[-1].lnk), (yyvsp[0].lnk), "type_qualifier_list type_qualifier skipped");
               }
#line 4787 "y.tab.c"
    break;

  case 222: /* type_qualifier_list_opt: %empty  */
#line 1486 "SDCC.y"
    {
      (yyval.lnk) = 0;
    }
#line 4795 "y.tab.c"
    break;

  case 223: /* type_qualifier_list_opt: type_qualifier_list  */
#line 1490 "SDCC.y"
    {
      (yyval.lnk) = (yyvsp[0].lnk);
    }
#line 4803 "y.tab.c"
    break;

  case 225: /* parameter_type_list: parameter_list ',' ELLIPSIS  */
#line 1497 "SDCC.y"
                                      { (yyvsp[-2].val)->vArgs = 1;}
#line 4809 "y.tab.c"
    break;

  case 227: /* parameter_list: parameter_list ',' parameter_declaration  */
#line 1503 "SDCC.y"
         {
            (yyvsp[0].val)->next = (yyvsp[-2].val);
            (yyval.val) = (yyvsp[0].val);
         }
#line 4818 "y.tab.c"
    break;

  case 228: /* parameter_declaration: declaration_specifiers declarator  */
#line 1511 "SDCC.y"
        {
          symbol *loop;

          if (IS_SPEC ((yyvsp[-1].lnk)) && !IS_VALID_PARAMETER_STORAGE_CLASS_SPEC ((yyvsp[-1].lnk)))
            {
              werror (E_STORAGE_CLASS_FOR_PARAMETER, (yyvsp[0].sym)->name);
            }
          pointerTypes ((yyvsp[0].sym)->type, (yyvsp[-1].lnk));
          if (IS_SPEC ((yyvsp[0].sym)->etype))
            SPEC_NEEDSPAR((yyvsp[0].sym)->etype) = 0;
          addDecl ((yyvsp[0].sym), 0, (yyvsp[-1].lnk));
          for (loop = (yyvsp[0].sym); loop; loop->_isparm = 1, loop = loop->next)
            ;
          (yyval.val) = symbolVal ((yyvsp[0].sym));
          ignoreTypedefType = 0;
        }
#line 4839 "y.tab.c"
    break;

  case 229: /* parameter_declaration: type_name  */
#line 1528 "SDCC.y"
        {
          (yyval.val) = newValue ();
          (yyval.val)->type = (yyvsp[0].lnk);
          (yyval.val)->etype = getSpec ((yyval.val)->type);
          ignoreTypedefType = 0;
         }
#line 4850 "y.tab.c"
    break;

  case 230: /* abstract_declarator: pointer  */
#line 1537 "SDCC.y"
             { (yyval.lnk) = reverseLink((yyvsp[0].lnk)); }
#line 4856 "y.tab.c"
    break;

  case 232: /* abstract_declarator: pointer direct_abstract_declarator  */
#line 1539 "SDCC.y"
                                          { (yyvsp[-1].lnk) = reverseLink((yyvsp[-1].lnk)); (yyvsp[0].lnk)->next = (yyvsp[-1].lnk); (yyval.lnk) = (yyvsp[0].lnk);
          if (IS_PTR((yyvsp[-1].lnk)) && IS_FUNC((yyvsp[0].lnk)))
            DCL_TYPE((yyvsp[-1].lnk)) = CPOINTER;
        }
#line 4865 "y.tab.c"
    break;

  case 233: /* direct_abstract_declarator: '(' abstract_declarator ')'  */
#line 1546 "SDCC.y"
                                    { (yyval.lnk) = (yyvsp[-1].lnk); }
#line 4871 "y.tab.c"
    break;

  case 236: /* direct_abstract_declarator_opt: %empty  */
#line 1552 "SDCC.y"
                                 { (yyval.lnk) = NULL; }
#line 4877 "y.tab.c"
    break;

  case 238: /* array_abstract_declarator: direct_abstract_declarator_opt '[' ']'  */
#line 1557 "SDCC.y"
                                              {
                                       (yyval.lnk) = newLink (DECLARATOR);
                                       DCL_TYPE((yyval.lnk)) = ARRAY;
                                       DCL_ELEM((yyval.lnk)) = 0;
                                       if((yyvsp[-2].lnk))
                                         (yyval.lnk)->next = (yyvsp[-2].lnk);
                                    }
#line 4889 "y.tab.c"
    break;

  case 239: /* array_abstract_declarator: direct_abstract_declarator_opt '[' constant_expr ']'  */
#line 1565 "SDCC.y"
                                    {
                                       value *val;
                                       (yyval.lnk) = newLink (DECLARATOR);
                                       DCL_TYPE((yyval.lnk)) = ARRAY;
                                       DCL_ELEM((yyval.lnk)) = (int) ulFromVal(val = constExprValue((yyvsp[-1].asts), true));
                                       if((yyvsp[-3].lnk))
                                         (yyval.lnk)->next = (yyvsp[-3].lnk);
                                    }
#line 4902 "y.tab.c"
    break;

  case 240: /* function_abstract_declarator: '(' ')'  */
#line 1576 "SDCC.y"
                                    { (yyval.lnk) = NULL;}
#line 4908 "y.tab.c"
    break;

  case 241: /* function_abstract_declarator: direct_abstract_declarator '(' ')'  */
#line 1577 "SDCC.y"
                                        {
     // $1 must be a pointer to a function
     sym_link *p=newLink(DECLARATOR);
     DCL_TYPE(p) = FUNCTION;
     if (!(yyvsp[-2].lnk)) {
       // ((void (code *) ()) 0) ()
       (yyvsp[-2].lnk)=newLink(DECLARATOR);
       DCL_TYPE((yyvsp[-2].lnk))=CPOINTER;
       (yyval.lnk) = (yyvsp[-2].lnk);
     }
     (yyvsp[-2].lnk)->next=p;
   }
#line 4925 "y.tab.c"
    break;

  case 242: /* function_abstract_declarator: '(' parameter_type_list ')'  */
#line 1589 "SDCC.y"
                                    { (yyval.lnk) = NULL;}
#line 4931 "y.tab.c"
    break;

  case 243: /* $@9: %empty  */
#line 1591 "SDCC.y"
        {
          NestLevel += LEVEL_UNIT;
          STACK_PUSH(blockNum, currBlockno);
          btree_add_child(currBlockno, ++blockNo);
          currBlockno = blockNo;
        }
#line 4942 "y.tab.c"
    break;

  case 244: /* function_abstract_declarator: direct_abstract_declarator '(' $@9 parameter_type_list ')'  */
#line 1598 "SDCC.y"
        {
          sym_link *p = newLink(DECLARATOR), *q;
          DCL_TYPE(p) = FUNCTION;

          FUNC_HASVARARGS(p) = IS_VARG((yyvsp[-1].val));
          FUNC_ARGS(p) = reverseVal((yyvsp[-1].val));

          /* nest level was incremented to take care of the parms  */
          NestLevel -= LEVEL_UNIT;
          currBlockno = STACK_POP(blockNum);
          if (!(yyvsp[-4].lnk))
            {
              /* ((void (code *) (void)) 0) () */
              (yyvsp[-4].lnk) = newLink(DECLARATOR);
              DCL_TYPE((yyvsp[-4].lnk)) = CPOINTER;
              (yyval.lnk) = (yyvsp[-4].lnk);
            }
          for (q = (yyvsp[-4].lnk); q && q->next; q = q->next);
          q->next = p;
        }
#line 4967 "y.tab.c"
    break;

  case 245: /* initializer: assignment_expr  */
#line 1621 "SDCC.y"
                                    { (yyval.ilist) = newiList(INIT_NODE,(yyvsp[0].asts)); }
#line 4973 "y.tab.c"
    break;

  case 246: /* initializer: '{' '}'  */
#line 1623 "SDCC.y"
     {
       if (!options.std_c2x)
         werror(W_EMPTY_INIT_C2X);
       // {} behaves mostly like {0}, so we emulate that, and use the isempty flag for the few cases where there is a difference.
       (yyval.ilist) = newiList(INIT_DEEP, revinit(newiList(INIT_NODE, newAst_VALUE(constIntVal("0")))));
       (yyval.ilist)->isempty = true;
     }
#line 4985 "y.tab.c"
    break;

  case 247: /* initializer: '{' initializer_list '}'  */
#line 1630 "SDCC.y"
                                    { (yyval.ilist) = newiList(INIT_DEEP,revinit((yyvsp[-1].ilist))); }
#line 4991 "y.tab.c"
    break;

  case 248: /* initializer: '{' initializer_list ',' '}'  */
#line 1631 "SDCC.y"
                                    { (yyval.ilist) = newiList(INIT_DEEP,revinit((yyvsp[-2].ilist))); }
#line 4997 "y.tab.c"
    break;

  case 249: /* initializer_list: designation_opt initializer  */
#line 1635 "SDCC.y"
                                    { (yyvsp[0].ilist)->designation = (yyvsp[-1].dsgn); (yyval.ilist) = (yyvsp[0].ilist); }
#line 5003 "y.tab.c"
    break;

  case 250: /* initializer_list: initializer_list ',' designation_opt initializer  */
#line 1637 "SDCC.y"
                                    {
                                       (yyvsp[0].ilist)->designation = (yyvsp[-1].dsgn);
                                       (yyvsp[0].ilist)->next = (yyvsp[-3].ilist);
                                       (yyval.ilist) = (yyvsp[0].ilist);
                                    }
#line 5013 "y.tab.c"
    break;

  case 251: /* designation_opt: %empty  */
#line 1645 "SDCC.y"
                                 { (yyval.dsgn) = NULL; }
#line 5019 "y.tab.c"
    break;

  case 253: /* designation: designator_list '='  */
#line 1650 "SDCC.y"
                                 { (yyval.dsgn) = revDesignation((yyvsp[-1].dsgn)); }
#line 5025 "y.tab.c"
    break;

  case 255: /* designator_list: designator_list designator  */
#line 1655 "SDCC.y"
                                 { (yyvsp[0].dsgn)->next = (yyvsp[-1].dsgn); (yyval.dsgn) = (yyvsp[0].dsgn); }
#line 5031 "y.tab.c"
    break;

  case 256: /* designator: '[' constant_expr ']'  */
#line 1660 "SDCC.y"
         {
            value *tval;
            int elemno;

            tval = constExprValue((yyvsp[-1].asts), true);
            /* if it is not a constant then Error  */
            if (!tval || (SPEC_SCLS(tval->etype) != S_LITERAL))
              {
                werror (E_CONST_EXPECTED);
                elemno = 0; /* arbitrary fixup */
              }
            else
              {
                if ((elemno = (int) ulFromVal(tval)) < 0)
                  {
                    werror (E_BAD_DESIGNATOR);
                    elemno = 0; /* arbitrary fixup */
                  }
              }
            (yyval.dsgn) = newDesignation(DESIGNATOR_ARRAY, &elemno);
         }
#line 5057 "y.tab.c"
    break;

  case 257: /* designator: '.' identifier  */
#line 1681 "SDCC.y"
                                 { (yyval.dsgn) = newDesignation(DESIGNATOR_STRUCT,(yyvsp[0].sym)); }
#line 5063 "y.tab.c"
    break;

  case 258: /* static_assert_declaration: STATIC_ASSERT '(' constant_expr ',' STRING_LITERAL ')' ';'  */
#line 1686 "SDCC.y"
                                    {
                                       value *val = constExprValue ((yyvsp[-4].asts), true);
                                       if (!val)
                                         werror (E_CONST_EXPECTED);
                                       else if (!ullFromVal(val))
                                         werror (W_STATIC_ASSERTION, (yyvsp[-2].yystr));
                                    }
#line 5075 "y.tab.c"
    break;

  case 259: /* static_assert_declaration: STATIC_ASSERT '(' constant_expr ')' ';'  */
#line 1694 "SDCC.y"
                                    {
                                       value *val = constExprValue ((yyvsp[-2].asts), true);
                                       if (!options.std_c2x)
                                         werror (E_STATIC_ASSERTION_C2X);
                                       if (!val)
                                         werror (E_CONST_EXPECTED);
                                       else if (!ullFromVal(val))
                                         werror (W_STATIC_ASSERTION_2);
                                    }
#line 5089 "y.tab.c"
    break;

  case 260: /* attribute_specifier_sequence: attribute_specifier_sequence attribute_specifier  */
#line 1707 "SDCC.y"
     {
       (yyval.attr) = (yyvsp[-1].attr);
       attribute *a;
       for (a = (yyval.attr); a->next; a = a->next);
       a->next = (yyvsp[0].attr);
     }
#line 5100 "y.tab.c"
    break;

  case 261: /* attribute_specifier_sequence: attribute_specifier  */
#line 1714 "SDCC.y"
     {
       (yyval.attr) = (yyvsp[0].attr);
     }
#line 5108 "y.tab.c"
    break;

  case 262: /* attribute_specifier_sequence_opt: %empty  */
#line 1721 "SDCC.y"
     {
       (yyval.attr) = 0;
     }
#line 5116 "y.tab.c"
    break;

  case 263: /* attribute_specifier_sequence_opt: attribute_specifier_sequence  */
#line 1725 "SDCC.y"
     {
       (yyval.attr) = (yyvsp[0].attr);
     }
#line 5124 "y.tab.c"
    break;

  case 264: /* attribute_specifier: ATTR_START attribute_list ']' ']'  */
#line 1732 "SDCC.y"
     {
       if (!options.std_c2x)
         werror(E_ATTRIBUTE_C2X);
       (yyval.attr) = (yyvsp[-2].attr);
     }
#line 5134 "y.tab.c"
    break;

  case 265: /* attribute_list: attribute_opt  */
#line 1741 "SDCC.y"
     {
       (yyval.attr) = (yyvsp[0].attr);
     }
#line 5142 "y.tab.c"
    break;

  case 266: /* attribute_list: attribute_list ',' attribute_opt  */
#line 1745 "SDCC.y"
     {
       (yyval.attr) = (yyvsp[-2].attr);
       if ((yyvsp[0].attr))
         {
           attribute *a;
           for (a = (yyval.attr); a->next; a = a->next);
           a->next = (yyvsp[0].attr);
         }
     }
#line 5156 "y.tab.c"
    break;

  case 267: /* attribute: attribute_token  */
#line 1758 "SDCC.y"
   {
     (yyval.attr) = newAttribute ((yyvsp[0].sym), 0);
   }
#line 5164 "y.tab.c"
    break;

  case 268: /* attribute: attribute_token attribute_argument_clause  */
#line 1762 "SDCC.y"
   {
     (yyval.attr) = newAttribute ((yyvsp[-1].sym), 0);
   }
#line 5172 "y.tab.c"
    break;

  case 269: /* attribute_opt: %empty  */
#line 1769 "SDCC.y"
     {
       (yyval.attr) = 0;
     }
#line 5180 "y.tab.c"
    break;

  case 270: /* attribute_opt: attribute  */
#line 1773 "SDCC.y"
     {
       (yyval.attr) = (yyvsp[0].attr);
     }
#line 5188 "y.tab.c"
    break;

  case 271: /* attribute_token: identifier  */
#line 1780 "SDCC.y"
     {
       (yyval.sym) = (yyvsp[0].sym);
       (yyval.sym)->next = 0;
       werror (W_UNKNOWN_ATTRIBUTE, (yyvsp[0].sym)->name);
     }
#line 5198 "y.tab.c"
    break;

  case 272: /* attribute_token: identifier TOK_SEP identifier  */
#line 1786 "SDCC.y"
     {
       (yyval.sym) = (yyvsp[-2].sym);
       (yyval.sym)->next = (yyvsp[0].sym);
       struct dbuf_s dbuf;
       dbuf_init (&dbuf, 64);
       dbuf_printf (&dbuf, "%s::%s", (yyvsp[-2].sym)->name, (yyvsp[0].sym)->name);
       werror (W_UNKNOWN_ATTRIBUTE, dbuf_c_str (&dbuf));
       dbuf_destroy (&dbuf);
     }
#line 5212 "y.tab.c"
    break;

  case 283: /* unlabeled_statement: attribute_specifier_sequence_opt primary_block  */
#line 1824 "SDCC.y"
     {
       (yyval.asts) = (yyvsp[0].asts);
     }
#line 5220 "y.tab.c"
    break;

  case 284: /* unlabeled_statement: attribute_specifier_sequence_opt jump_statement  */
#line 1828 "SDCC.y"
     {
       (yyval.asts) = (yyvsp[0].asts);
     }
#line 5228 "y.tab.c"
    break;

  case 287: /* primary_block: compound_statement  */
#line 1837 "SDCC.y"
     {
       (yyval.asts) = (yyvsp[0].asts);
     }
#line 5236 "y.tab.c"
    break;

  case 288: /* primary_block: selection_statement  */
#line 1841 "SDCC.y"
     {
       (yyval.asts) = (yyvsp[0].asts);
     }
#line 5244 "y.tab.c"
    break;

  case 289: /* primary_block: iteration_statement  */
#line 1845 "SDCC.y"
     {
       (yyval.asts) = (yyvsp[0].asts);
     }
#line 5252 "y.tab.c"
    break;

  case 290: /* secondary_block: statement  */
#line 1852 "SDCC.y"
     {
       (yyval.asts) = (yyvsp[0].asts);
     }
#line 5260 "y.tab.c"
    break;

  case 291: /* labeled_statement: label statement  */
#line 1859 "SDCC.y"
     {
       if ((yyvsp[-1].asts))
         {
           (yyval.asts) = (yyvsp[-1].asts);
           (yyvsp[-1].asts)->right = (yyvsp[0].asts);
         }
       else
         (yyval.asts) = newNode (BLOCK, NULL, NULL);
     }
#line 5274 "y.tab.c"
    break;

  case 292: /* label: identifier ':'  */
#line 1872 "SDCC.y"
     {
       (yyval.asts) = createLabel((yyvsp[-1].sym),NULL);
       (yyvsp[-1].sym)->isitmp = 0;
     }
#line 5283 "y.tab.c"
    break;

  case 293: /* label: attribute_specifier_sequence_opt CASE constant_expr ':'  */
#line 1877 "SDCC.y"
     {
       if (STACK_EMPTY(swStk))
         (yyval.asts) = createCase(NULL,(yyvsp[-1].asts),NULL);
       else
         (yyval.asts) = createCase(STACK_PEEK(swStk),(yyvsp[-1].asts),NULL);
     }
#line 5294 "y.tab.c"
    break;

  case 294: /* @10: %empty  */
#line 1883 "SDCC.y"
                                              { (yyval.asts) = newNode(DEFAULT,NULL,NULL); }
#line 5300 "y.tab.c"
    break;

  case 295: /* label: attribute_specifier_sequence_opt DEFAULT @10 ':'  */
#line 1884 "SDCC.y"
     {
       if (STACK_EMPTY(swStk))
         (yyval.asts) = createDefault(NULL,(yyvsp[-1].asts),NULL);
       else
         (yyval.asts) = createDefault(STACK_PEEK(swStk),(yyvsp[-1].asts),NULL);
     }
#line 5311 "y.tab.c"
    break;

  case 296: /* start_block: '{'  */
#line 1894 "SDCC.y"
        {
          NestLevel += LEVEL_UNIT;
          STACK_PUSH(blockNum, currBlockno);
          btree_add_child(currBlockno, ++blockNo);
          currBlockno = blockNo;
          ignoreTypedefType = 0;
        }
#line 5323 "y.tab.c"
    break;

  case 297: /* end_block: '}'  */
#line 1905 "SDCC.y"
        {
          leaveBlockScope (currBlockno);
          NestLevel -= LEVEL_UNIT;
          currBlockno = STACK_POP(blockNum);
        }
#line 5333 "y.tab.c"
    break;

  case 298: /* compound_statement: start_block end_block  */
#line 1913 "SDCC.y"
                                              { (yyval.asts) = createBlock(NULL, NULL); }
#line 5339 "y.tab.c"
    break;

  case 299: /* compound_statement: start_block block_item_list end_block  */
#line 1915 "SDCC.y"
     {
       (yyval.asts) = (yyvsp[-1].asts);
       cleanUpLevel(StructTab, NestLevel + LEVEL_UNIT);
     }
#line 5348 "y.tab.c"
    break;

  case 300: /* block_item_list: statements_and_implicit  */
#line 1922 "SDCC.y"
                                              { (yyval.asts) = createBlock(NULL, (yyvsp[0].asts)); }
#line 5354 "y.tab.c"
    break;

  case 301: /* block_item_list: declaration_list  */
#line 1923 "SDCC.y"
                                              { (yyval.asts) = createBlock((yyvsp[0].sym), NULL); }
#line 5360 "y.tab.c"
    break;

  case 302: /* block_item_list: declaration_list statements_and_implicit  */
#line 1924 "SDCC.y"
                                              { (yyval.asts) = createBlock((yyvsp[-1].sym), (yyvsp[0].asts)); }
#line 5366 "y.tab.c"
    break;

  case 304: /* expression_statement: attribute_specifier_sequence expression ';'  */
#line 1929 "SDCC.y"
                                                           { (yyval.asts) = (yyvsp[-1].asts); seqPointNo++;}
#line 5372 "y.tab.c"
    break;

  case 305: /* else_statement: ELSE statement  */
#line 1933 "SDCC.y"
                        { (yyval.asts) = (yyvsp[0].asts); }
#line 5378 "y.tab.c"
    break;

  case 306: /* else_statement: %empty  */
#line 1934 "SDCC.y"
                        { (yyval.asts) = NULL; }
#line 5384 "y.tab.c"
    break;

  case 307: /* $@11: %empty  */
#line 1939 "SDCC.y"
                           { seqPointNo++;}
#line 5390 "y.tab.c"
    break;

  case 308: /* selection_statement: IF '(' expression ')' $@11 statement else_statement  */
#line 1940 "SDCC.y"
                           {
                              noLineno++;
                              (yyval.asts) = createIf ((yyvsp[-4].asts), (yyvsp[-1].asts), (yyvsp[0].asts) );
                              (yyval.asts)->lineno = (yyvsp[-4].asts)->lineno;
                              (yyval.asts)->filename = (yyvsp[-4].asts)->filename;
                              noLineno--;
                           }
#line 5402 "y.tab.c"
    break;

  case 309: /* @12: %empty  */
#line 1947 "SDCC.y"
                                 {
                              ast *ex;
                              static   int swLabel = 0;

                              seqPointNo++;
                              /* create a node for expression  */
                              ex = newNode(SWITCH,(yyvsp[-1].asts),NULL);
                              STACK_PUSH(swStk,ex);   /* save it in the stack */
                              ex->values.switchVals.swNum = swLabel;

                              /* now create the label */
                              SNPRINTF(lbuff, sizeof(lbuff),
                                       "_swBrk_%d",swLabel++);
                              (yyval.sym)  =  newSymbol(lbuff,NestLevel);
                              /* put label in the break stack  */
                              STACK_PUSH(breakStack,(yyval.sym));
                           }
#line 5424 "y.tab.c"
    break;

  case 310: /* selection_statement: SWITCH '(' expression ')' @12 secondary_block  */
#line 1964 "SDCC.y"
                           {
                              /* get back the switch form the stack  */
                              (yyval.asts) = STACK_POP(swStk);
                              (yyval.asts)->right = newNode (NULLOP,(yyvsp[0].asts),createLabel((yyvsp[-1].sym),NULL));
                              STACK_POP(breakStack);
                           }
#line 5435 "y.tab.c"
    break;

  case 311: /* $@13: %empty  */
#line 1973 "SDCC.y"
                              { seqPointNo++;}
#line 5441 "y.tab.c"
    break;

  case 312: /* iteration_statement: while '(' expression ')' $@13 secondary_block  */
#line 1974 "SDCC.y"
                         {
                           noLineno++;
                           (yyval.asts) = createWhile ( (yyvsp[-5].sym), STACK_POP(continueStack),
                                              STACK_POP(breakStack), (yyvsp[-3].asts), (yyvsp[0].asts) );
                           (yyval.asts)->lineno = (yyvsp[-5].sym)->lineDef;
                           (yyval.asts)->filename = (yyvsp[-5].sym)->fileDef;
                           noLineno--;
                         }
#line 5454 "y.tab.c"
    break;

  case 313: /* iteration_statement: do secondary_block WHILE '(' expression ')' ';'  */
#line 1983 "SDCC.y"
                        {
                          seqPointNo++;
                          noLineno++;
                          (yyval.asts) = createDo ( (yyvsp[-6].sym) , STACK_POP(continueStack),
                                          STACK_POP(breakStack), (yyvsp[-2].asts), (yyvsp[-5].asts));
                          (yyval.asts)->lineno = (yyvsp[-6].sym)->lineDef;
                          (yyval.asts)->filename = (yyvsp[-6].sym)->fileDef;
                          noLineno--;
                        }
#line 5468 "y.tab.c"
    break;

  case 314: /* iteration_statement: for '(' expression_opt ';' expression_opt ';' expression_opt ')' secondary_block  */
#line 1993 "SDCC.y"
                        {
                          noLineno++;

                          (yyval.asts) = newNode(FOR,(yyvsp[0].asts),NULL);
                          AST_FOR((yyval.asts),trueLabel) = (yyvsp[-8].sym);
                          AST_FOR((yyval.asts),continueLabel) =  STACK_POP(continueStack);
                          AST_FOR((yyval.asts),falseLabel) = STACK_POP(breakStack);
                          AST_FOR((yyval.asts),condLabel)  = STACK_POP(forStack);
                          AST_FOR((yyval.asts),initExpr)   = (yyvsp[-6].asts);
                          AST_FOR((yyval.asts),condExpr)   = (yyvsp[-4].asts);
                          AST_FOR((yyval.asts),loopExpr)   = (yyvsp[-2].asts);
                          
                          /* This continue label is not at the correct location, */
                          /* but we need to create it now for proper binding. The */
                          /* code that handles the FOR node will move it to the */
                          /* proper location inside the for loop. */
                          if (AST_FOR((yyval.asts),continueLabel)->isref)
                            (yyval.asts)->right = createLabel(AST_FOR((yyval.asts),continueLabel),NULL);
                          (yyval.asts) = newNode(NULLOP,(yyval.asts),createLabel(AST_FOR((yyval.asts),falseLabel),NULL));
                          noLineno--;

                          NestLevel -= LEVEL_UNIT;
                          currBlockno = STACK_POP(blockNum);
                        }
#line 5497 "y.tab.c"
    break;

  case 315: /* $@14: %empty  */
#line 2018 "SDCC.y"
                        {
                          if (!options.std_c99)
                            werror (E_FOR_INITAL_DECLARATION_C99);

                          if ( (yyvsp[-4].sym) && IS_TYPEDEF((yyvsp[-4].sym)->etype))
                            allocVariables ((yyvsp[-4].sym));
                          ignoreTypedefType = 0;
                          addSymChain(&(yyvsp[-4].sym));
                        }
#line 5511 "y.tab.c"
    break;

  case 316: /* iteration_statement: for '(' declaration expression_opt ';' expression_opt ')' $@14 secondary_block  */
#line 2028 "SDCC.y"
                        {

                          noLineno++;

                          (yyval.asts) = newNode(FOR,(yyvsp[0].asts),NULL);
                          AST_FOR((yyval.asts),trueLabel) = (yyvsp[-8].sym);
                          AST_FOR((yyval.asts),continueLabel) =  STACK_POP(continueStack);
                          AST_FOR((yyval.asts),falseLabel) = STACK_POP(breakStack);
                          AST_FOR((yyval.asts),condLabel)  = STACK_POP(forStack);
                          AST_FOR((yyval.asts),initExpr)   = 0;
                          AST_FOR((yyval.asts),condExpr)   = (yyvsp[-5].asts);
                          AST_FOR((yyval.asts),loopExpr)   = (yyvsp[-3].asts);

                          /* This continue label is not at the correct location, */
                          /* but we need to create it now for proper binding. The */
                          /* code that handles the FOR node will move it to the */
                          /* proper location inside the for loop. */
                          if (AST_FOR((yyval.asts),continueLabel)->isref)
                            (yyval.asts)->right = createLabel(AST_FOR((yyval.asts),continueLabel),NULL);
                          (yyval.asts) = createBlock((yyvsp[-6].sym), newNode(NULLOP,(yyval.asts),createLabel(AST_FOR((yyval.asts),falseLabel),NULL)));
                          cleanUpLevel(StructTab, NestLevel + LEVEL_UNIT);
                          noLineno--;

                          leaveBlockScope (currBlockno);
                          NestLevel -= LEVEL_UNIT;
                          currBlockno = STACK_POP(blockNum);
                        }
#line 5543 "y.tab.c"
    break;

  case 317: /* jump_statement: GOTO identifier ';'  */
#line 2058 "SDCC.y"
                           {
                              if (inCriticalBlock) {
                                werror(E_INVALID_CRITICAL);
                                (yyval.asts) = NULL;
                              } else {
                                (yyvsp[-1].sym)->islbl = 1;
                                (yyval.asts) = newAst_VALUE(symbolVal((yyvsp[-1].sym)));
                                (yyval.asts) = newNode(GOTO,(yyval.asts),NULL);
                              }
                           }
#line 5558 "y.tab.c"
    break;

  case 318: /* jump_statement: CONTINUE ';'  */
#line 2068 "SDCC.y"
                           {
       /* make sure continue is in context */
       if (STACK_EMPTY(continueStack) || STACK_PEEK(continueStack) == NULL) {
           werror(E_BREAK_CONTEXT);
           (yyval.asts) = NULL;
       }
       else {
           (yyval.asts) = newAst_VALUE(symbolVal(STACK_PEEK(continueStack)));
           (yyval.asts) = newNode(GOTO,(yyval.asts),NULL);
           /* mark the continue label as referenced */
           STACK_PEEK(continueStack)->isref = 1;
       }
   }
#line 5576 "y.tab.c"
    break;

  case 319: /* jump_statement: BREAK ';'  */
#line 2081 "SDCC.y"
                           {
       if (STACK_EMPTY(breakStack) || STACK_PEEK(breakStack) == NULL) {
           werror(E_BREAK_CONTEXT);
           (yyval.asts) = NULL;
       } else {
           (yyval.asts) = newAst_VALUE(symbolVal(STACK_PEEK(breakStack)));
           (yyval.asts) = newNode(GOTO,(yyval.asts),NULL);
           STACK_PEEK(breakStack)->isref = 1;
       }
   }
#line 5591 "y.tab.c"
    break;

  case 320: /* jump_statement: RETURN ';'  */
#line 2091 "SDCC.y"
                           {
       seqPointNo++;
       if (inCriticalBlock) {
           werror(E_INVALID_CRITICAL);
           (yyval.asts) = NULL;
       } else {
           (yyval.asts) = newNode(RETURN,NULL,NULL);
       }
   }
#line 5605 "y.tab.c"
    break;

  case 321: /* jump_statement: RETURN expression ';'  */
#line 2100 "SDCC.y"
                                 {
       seqPointNo++;
       if (inCriticalBlock) {
           werror(E_INVALID_CRITICAL);
           (yyval.asts) = NULL;
       } else {
           (yyval.asts) = newNode(RETURN,NULL,(yyvsp[-1].asts));
       }
   }
#line 5619 "y.tab.c"
    break;

  case 324: /* external_declaration: function_definition  */
#line 2120 "SDCC.y"
        {
          // blockNo = 0;
        }
#line 5627 "y.tab.c"
    break;

  case 325: /* external_declaration: declaration  */
#line 2124 "SDCC.y"
        {
          ignoreTypedefType = 0;
          if ((yyvsp[0].sym) && (yyvsp[0].sym)->type && IS_FUNC((yyvsp[0].sym)->type))
            {
              /* The only legal storage classes for
               * a function prototype (declaration)
               * are extern and static. extern is the
               * default. Thus, if this function isn't
               * explicitly marked static, mark it
               * extern.
               */
              if ((yyvsp[0].sym)->etype && IS_SPEC((yyvsp[0].sym)->etype) && !SPEC_STAT((yyvsp[0].sym)->etype))
                {
                  SPEC_EXTR((yyvsp[0].sym)->etype) = 1;
                }
            }
          addSymChain (&(yyvsp[0].sym));
          allocVariables ((yyvsp[0].sym));
          cleanUpLevel (SymbolTab, 1);
        }
#line 5652 "y.tab.c"
    break;

  case 326: /* external_declaration: addressmod  */
#line 2145 "SDCC.y"
        {
          /* These empty braces here are apparently required by some version of GNU bison on MS Windows. See bug #2858. */
        }
#line 5660 "y.tab.c"
    break;

  case 327: /* $@15: %empty  */
#line 2152 "SDCC.y"
        {   /* function return type not specified, which is allowed in C90 (and means int), but disallowed in C99 and later */
            werror (W_RETURN_TYPE_OMITTED_INT);
            addDecl((yyvsp[0].sym),0,newIntLink());
            (yyvsp[0].sym) = createFunctionDecl((yyvsp[0].sym));
            if ((yyvsp[0].sym))
                {
                    if (FUNC_ISCRITICAL ((yyvsp[0].sym)->type))
                        inCriticalFunction = 1;
                    strncpy (function_name, (yyvsp[0].sym)->name, sizeof (function_name) - 4);
                    memset (function_name + sizeof (function_name) - 4, 0x00, 4);
                }
        }
#line 5677 "y.tab.c"
    break;

  case 328: /* function_definition: declarator $@15 function_body  */
#line 2165 "SDCC.y"
        {
            (yyval.asts) = createFunction((yyvsp[-2].sym),(yyvsp[0].asts));
            if ((yyvsp[-2].sym) && FUNC_ISCRITICAL ((yyvsp[-2].sym)->type))
                inCriticalFunction = 0;
        }
#line 5687 "y.tab.c"
    break;

  case 329: /* $@16: %empty  */
#line 2171 "SDCC.y"
        {
            sym_link *p = copyLinkChain((yyvsp[-1].lnk));
            pointerTypes((yyvsp[0].sym)->type,p);
            addDecl((yyvsp[0].sym),0,p);
            (yyvsp[0].sym) = createFunctionDecl((yyvsp[0].sym));
            if ((yyvsp[0].sym))
                {
                	if (!strcmp ((yyvsp[0].sym)->name, "_sdcc_external_startup")) // The rename (and semantics change happened) in SDCC 4.2.10. Keep this warning for two major releases afterwards.
                		werror (W__SDCC_EXTERNAL_STARTUP_DEF);
                    if (FUNC_ISCRITICAL ((yyvsp[0].sym)->type))
                        inCriticalFunction = 1;
                    // warn for loss of calling convention for inlined functions.
                    if (FUNC_ISINLINE ((yyvsp[0].sym)->type) &&
                        ( FUNC_ISZ88DK_CALLEE ((yyvsp[0].sym)->type) || FUNC_ISZ88DK_FASTCALL ((yyvsp[0].sym)->type) ||
                          FUNC_BANKED ((yyvsp[0].sym)->type)         || FUNC_REGBANK ((yyvsp[0].sym)->type)          ||
                          FUNC_ISOVERLAY ((yyvsp[0].sym)->type)      || FUNC_ISISR ((yyvsp[0].sym)->type) ))
                        werror (W_INLINE_FUNCATTR, (yyvsp[0].sym)->name);
                    strncpy (function_name, (yyvsp[0].sym)->name, sizeof (function_name) - 4);
                    memset (function_name + sizeof (function_name) - 4, 0x00, 4);
                }
        }
#line 5713 "y.tab.c"
    break;

  case 330: /* function_definition: declaration_specifiers declarator $@16 function_body  */
#line 2193 "SDCC.y"
        {
            (yyval.asts) = createFunction((yyvsp[-2].sym),(yyvsp[0].asts));
            if ((yyvsp[-2].sym) && FUNC_ISCRITICAL ((yyvsp[-2].sym)->type))
                inCriticalFunction = 0;
        }
#line 5723 "y.tab.c"
    break;

  case 332: /* function_body: kr_declaration_list compound_statement  */
#line 2203 "SDCC.y"
     {
       werror (E_OLD_STYLE, ((yyvsp[-1].sym) ? (yyvsp[-1].sym)->name: ""));
       exit (1);
     }
#line 5732 "y.tab.c"
    break;

  case 333: /* file: %empty  */
#line 2213 "SDCC.y"
     {
       werror(W_EMPTY_TRANSLATION_UNIT);
     }
#line 5740 "y.tab.c"
    break;

  case 336: /* function_attributes: function_attributes function_attribute  */
#line 2224 "SDCC.y"
                                            { (yyval.lnk) = mergeSpec((yyvsp[-1].lnk),(yyvsp[0].lnk),"function_attribute"); }
#line 5746 "y.tab.c"
    break;

  case 337: /* function_attribute: USING '(' constant_expr ')'  */
#line 2228 "SDCC.y"
                                  {
                        (yyval.lnk) = newLink(SPECIFIER);
                        FUNC_REGBANK((yyval.lnk)) = (int) ulFromVal(constExprValue((yyvsp[-1].asts), true));
                     }
#line 5755 "y.tab.c"
    break;

  case 338: /* function_attribute: REENTRANT  */
#line 2232 "SDCC.y"
                     {  (yyval.lnk) = newLink (SPECIFIER);
                        FUNC_ISREENT((yyval.lnk))=1;
                     }
#line 5763 "y.tab.c"
    break;

  case 339: /* function_attribute: CRITICAL  */
#line 2235 "SDCC.y"
                     {  (yyval.lnk) = newLink (SPECIFIER);
                        FUNC_ISCRITICAL((yyval.lnk)) = 1;
                     }
#line 5771 "y.tab.c"
    break;

  case 340: /* function_attribute: NAKED  */
#line 2238 "SDCC.y"
                     {  (yyval.lnk) = newLink (SPECIFIER);
                        FUNC_ISNAKED((yyval.lnk))=1;
                     }
#line 5779 "y.tab.c"
    break;

  case 341: /* function_attribute: JAVANATIVE  */
#line 2241 "SDCC.y"
                     {  (yyval.lnk) = newLink (SPECIFIER);
                        FUNC_ISJAVANATIVE((yyval.lnk))=1;
                     }
#line 5787 "y.tab.c"
    break;

  case 342: /* function_attribute: OVERLAY  */
#line 2244 "SDCC.y"
                     {  (yyval.lnk) = newLink (SPECIFIER);
                        FUNC_ISOVERLAY((yyval.lnk))=1;
                     }
#line 5795 "y.tab.c"
    break;

  case 343: /* function_attribute: NONBANKED  */
#line 2247 "SDCC.y"
                     {(yyval.lnk) = newLink (SPECIFIER);
                        FUNC_NONBANKED((yyval.lnk)) = 1;
                        if (FUNC_BANKED((yyval.lnk))) {
                            werror(W_BANKED_WITH_NONBANKED);
                        }
                     }
#line 5806 "y.tab.c"
    break;

  case 344: /* function_attribute: SHADOWREGS  */
#line 2253 "SDCC.y"
                     {(yyval.lnk) = newLink (SPECIFIER);
                        FUNC_ISSHADOWREGS((yyval.lnk)) = 1;
                     }
#line 5814 "y.tab.c"
    break;

  case 345: /* function_attribute: SD_WPARAM  */
#line 2256 "SDCC.y"
                     {(yyval.lnk) = newLink (SPECIFIER);
                        FUNC_ISWPARAM((yyval.lnk)) = 1;
                     }
#line 5822 "y.tab.c"
    break;

  case 346: /* function_attribute: BANKED  */
#line 2259 "SDCC.y"
                     {(yyval.lnk) = newLink (SPECIFIER);
                        FUNC_BANKED((yyval.lnk)) = 1;
                        if (FUNC_NONBANKED((yyval.lnk))) {
                            werror(W_BANKED_WITH_NONBANKED);
                        }
                     }
#line 5833 "y.tab.c"
    break;

  case 347: /* function_attribute: Interrupt_storage  */
#line 2266 "SDCC.y"
                     {
                        (yyval.lnk) = newLink (SPECIFIER);
                        FUNC_INTNO((yyval.lnk)) = (yyvsp[0].yyint);
                        FUNC_ISISR((yyval.lnk)) = 1;
                     }
#line 5843 "y.tab.c"
    break;

  case 348: /* function_attribute: TRAP  */
#line 2272 "SDCC.y"
                     {
                        (yyval.lnk) = newLink (SPECIFIER);
                        FUNC_INTNO((yyval.lnk)) = INTNO_TRAP;
                        FUNC_ISISR((yyval.lnk)) = 1;
                     }
#line 5853 "y.tab.c"
    break;

  case 349: /* function_attribute: SMALLC  */
#line 2277 "SDCC.y"
                     {  (yyval.lnk) = newLink (SPECIFIER);
                        FUNC_ISSMALLC((yyval.lnk)) = 1;
                     }
#line 5861 "y.tab.c"
    break;

  case 350: /* function_attribute: RAISONANCE  */
#line 2280 "SDCC.y"
                     {  (yyval.lnk) = newLink (SPECIFIER);
                        FUNC_ISRAISONANCE((yyval.lnk)) = 1;
                     }
#line 5869 "y.tab.c"
    break;

  case 351: /* function_attribute: IAR  */
#line 2283 "SDCC.y"
                     {  (yyval.lnk) = newLink (SPECIFIER);
                        FUNC_ISIAR((yyval.lnk)) = 1;
                     }
#line 5877 "y.tab.c"
    break;

  case 352: /* function_attribute: COSMIC  */
#line 2286 "SDCC.y"
                     {  (yyval.lnk) = newLink (SPECIFIER);
                        FUNC_ISCOSMIC((yyval.lnk)) = 1;
                     }
#line 5885 "y.tab.c"
    break;

  case 353: /* function_attribute: SDCCCALL '(' constant_expr ')'  */
#line 2290 "SDCC.y"
                     {  (yyval.lnk) = newLink (SPECIFIER);
                        FUNC_SDCCCALL((yyval.lnk)) = ulFromVal(constExprValue ((yyvsp[-1].asts), true));
                     }
#line 5893 "y.tab.c"
    break;

  case 354: /* function_attribute: Z88DK_FASTCALL  */
#line 2293 "SDCC.y"
                     {  (yyval.lnk) = newLink (SPECIFIER);
                        FUNC_ISZ88DK_FASTCALL((yyval.lnk)) = 1;
                     }
#line 5901 "y.tab.c"
    break;

  case 355: /* function_attribute: Z88DK_CALLEE  */
#line 2296 "SDCC.y"
                     {  (yyval.lnk) = newLink (SPECIFIER);
                        FUNC_ISZ88DK_CALLEE((yyval.lnk)) = 1;
                     }
#line 5909 "y.tab.c"
    break;

  case 356: /* function_attribute: Z88DK_PARAMS_OFFSET '(' constant_expr ')'  */
#line 2300 "SDCC.y"
                     {
                        value *offset_v = constExprValue ((yyvsp[-1].asts), true);
                        int    offset = 0;
                        (yyval.lnk) = newLink(SPECIFIER);
                        if  ( offset_v ) 
                          offset = ulFromVal(offset_v);
                        (yyval.lnk)->funcAttrs.z88dk_params_offset = offset;
                     }
#line 5922 "y.tab.c"
    break;

  case 357: /* function_attribute: Z88DK_SHORTCALL '(' constant_expr ',' constant_expr ')'  */
#line 2309 "SDCC.y"
                     {
                        value *rst_v = constExprValue ((yyvsp[-3].asts), true);
                        value *value_v = constExprValue ((yyvsp[-1].asts), true);
                        int rst = -1, value = -1;
                        (yyval.lnk) = newLink(SPECIFIER);

                        if  ( rst_v ) 
                          rst = ulFromVal(rst_v);
                        if  ( value_v ) 
                          value = ulFromVal(value_v);
          
                        if ( rst < 0 || rst > 56 || ( rst % 8 ) )
                          {
                            werror(E_SHORTCALL_INVALID_VALUE, "rst", rst);
                          }
                        if ( value < 0 || value > 0xfff )
                          {
                            werror(E_SHORTCALL_INVALID_VALUE, "value", value);
                          }
                        (yyval.lnk)->funcAttrs.z88dk_shortcall_rst = rst;
                        (yyval.lnk)->funcAttrs.z88dk_shortcall_val = value;
                        FUNC_ISZ88DK_SHORTCALL((yyval.lnk)) = 1;
                     }
#line 5950 "y.tab.c"
    break;

  case 358: /* function_attribute: PRESERVES_REGS '(' identifier_list ')'  */
#line 2333 "SDCC.y"
                     {
                        const struct symbol *regsym;
                        (yyval.lnk) = newLink (SPECIFIER);

                        for(regsym = (yyvsp[-1].sym); regsym; regsym = regsym->next)
                          {
                            int regnum;

                            if (!port->getRegByName || ((regnum = port->getRegByName(regsym->name)) < 0))
                              werror (W_UNKNOWN_REG, regsym->name);
                            else
                              (yyval.lnk)->funcAttrs.preserved_regs[regnum] = true;
                          }
                     }
#line 5969 "y.tab.c"
    break;

  case 359: /* offsetof_member_designator: identifier  */
#line 2350 "SDCC.y"
                     { (yyval.asts) = newAst_VALUE (symbolVal ((yyvsp[0].sym))); }
#line 5975 "y.tab.c"
    break;

  case 360: /* $@17: %empty  */
#line 2351 "SDCC.y"
                                    { ignoreTypedefType = 1; }
#line 5981 "y.tab.c"
    break;

  case 361: /* offsetof_member_designator: offsetof_member_designator '.' $@17 identifier  */
#line 2352 "SDCC.y"
                     {
                       ignoreTypedefType = 0;
                       (yyvsp[0].sym) = newSymbol ((yyvsp[0].sym)->name, NestLevel);
                       (yyvsp[0].sym)->implicit = 1;
                       (yyval.asts) = newNode ('.', (yyvsp[-3].asts), newAst_VALUE (symbolVal ((yyvsp[0].sym))));
                     }
#line 5992 "y.tab.c"
    break;

  case 362: /* offsetof_member_designator: offsetof_member_designator '[' expression ']'  */
#line 2359 "SDCC.y"
                     {
                       (yyval.asts) = newNode ('[', (yyvsp[-3].asts), (yyvsp[-1].asts));
                     }
#line 6000 "y.tab.c"
    break;

  case 363: /* string_literal_val: FUNC  */
#line 2366 "SDCC.y"
       {
         // essentially do $$ = newAst_VALUE (strVal("\"$function_name\""));

         value* val = newValue ();
         { // BUG: duplicate from strVal
           val->type = newLink (DECLARATOR);
           DCL_TYPE (val->type) = ARRAY;
           val->type->next = val->etype = newLink (SPECIFIER);
           SPEC_SCLS (val->etype) = S_LITERAL;
           SPEC_CONST (val->etype) = 1;
           SPEC_NOUN (val->etype) = V_CHAR;
           SPEC_USIGN (val->etype) = !options.signed_char;
           val->etype->select.s.b_implicit_sign = true;
         }

         int ll = 1 + strlen(function_name);
         char* s = (char*) Safe_alloc(ll*sizeof(char));
         if(ll > 1){
           s = strncpy(s, function_name, ll);
         }else{
           *s = 0;
         }
         SPEC_CVAL (val->etype).v_char = s;
         DCL_ELEM (val->type) = ll;
         (yyval.asts) = newAst_VALUE ( val );
       }
#line 6031 "y.tab.c"
    break;

  case 364: /* string_literal_val: STRING_LITERAL  */
#line 2392 "SDCC.y"
                    { (yyval.asts) = newAst_VALUE (strVal ((yyvsp[0].yystr))); }
#line 6037 "y.tab.c"
    break;

  case 365: /* Interrupt_storage: INTERRUPT  */
#line 2396 "SDCC.y"
               { (yyval.yyint) = INTNO_UNSPEC; }
#line 6043 "y.tab.c"
    break;

  case 366: /* Interrupt_storage: INTERRUPT '(' constant_expr ')'  */
#line 2398 "SDCC.y"
        { 
          value *val = constExprValue((yyvsp[-1].asts), true);
          int intno = (int) ulFromVal(val);
          if (val && (intno >= 0) && (intno <= INTNO_MAX))
            (yyval.yyint) = intno;
          else
            {
              werror(E_INT_BAD_INTNO, intno);
              (yyval.yyint) = INTNO_UNSPEC;
            }
        }
#line 6059 "y.tab.c"
    break;

  case 367: /* sfr_reg_bit: SBIT  */
#line 2412 "SDCC.y"
            {
               (yyval.lnk) = newLink(SPECIFIER);
               SPEC_NOUN((yyval.lnk)) = V_SBIT;
               SPEC_SCLS((yyval.lnk)) = S_SBIT;
               SPEC_BLEN((yyval.lnk)) = 1;
               SPEC_BSTR((yyval.lnk)) = 0;
               ignoreTypedefType = 1;
            }
#line 6072 "y.tab.c"
    break;

  case 369: /* sfr_attributes: SFR  */
#line 2424 "SDCC.y"
            {
               (yyval.lnk) = newLink(SPECIFIER);
               FUNC_REGBANK((yyval.lnk)) = 0;
               SPEC_NOUN((yyval.lnk))    = V_CHAR;
               SPEC_SCLS((yyval.lnk))    = S_SFR;
               SPEC_USIGN((yyval.lnk))   = 1;
               ignoreTypedefType = 1;
            }
#line 6085 "y.tab.c"
    break;

  case 370: /* sfr_attributes: SFR BANKED  */
#line 2432 "SDCC.y"
                {
               (yyval.lnk) = newLink(SPECIFIER);
               FUNC_REGBANK((yyval.lnk)) = 1;
               SPEC_NOUN((yyval.lnk))    = V_CHAR;
               SPEC_SCLS((yyval.lnk))    = S_SFR;
               SPEC_USIGN((yyval.lnk))   = 1;
               ignoreTypedefType = 1;
            }
#line 6098 "y.tab.c"
    break;

  case 371: /* sfr_attributes: SFR16  */
#line 2443 "SDCC.y"
            {
               (yyval.lnk) = newLink(SPECIFIER);
               FUNC_REGBANK((yyval.lnk)) = 0;
               SPEC_NOUN((yyval.lnk))    = V_INT;
               SPEC_SCLS((yyval.lnk))    = S_SFR;
               SPEC_USIGN((yyval.lnk))   = 1;
               ignoreTypedefType = 1;
            }
#line 6111 "y.tab.c"
    break;

  case 372: /* sfr_attributes: SFR32  */
#line 2454 "SDCC.y"
            {
               (yyval.lnk) = newLink(SPECIFIER);
               FUNC_REGBANK((yyval.lnk)) = 0;
               SPEC_NOUN((yyval.lnk))    = V_INT;
               SPEC_SCLS((yyval.lnk))    = S_SFR;
               SPEC_LONG((yyval.lnk))    = 1;
               SPEC_USIGN((yyval.lnk))   = 1;
               ignoreTypedefType = 1;
            }
#line 6125 "y.tab.c"
    break;

  case 374: /* opt_stag: %empty  */
#line 2467 "SDCC.y"
        {  /* synthesize a name add to structtable */
          ignoreTypedefType = 0;
          (yyval.sdef) = newStruct(genSymName(NestLevel));
          (yyval.sdef)->level = NestLevel;
          (yyval.sdef)->block = currBlockno;
          (yyval.sdef)->tagsym = NULL;
          //addSym (StructTab, $$, $$->tag, $$->level, currBlockno, 0);
        }
#line 6138 "y.tab.c"
    break;

  case 375: /* stag: identifier  */
#line 2479 "SDCC.y"
        {  /* add name to structure table */
          ignoreTypedefType = 0;
          (yyval.sdef) = newStruct((yyvsp[0].sym)->name);
          (yyval.sdef)->level = NestLevel;
          (yyval.sdef)->block = currBlockno;
          (yyval.sdef)->tagsym = (yyvsp[0].sym);
          //$$ = findSymWithBlock (StructTab, $1, currBlockno);
          //if (! $$ )
          //  {
          //    $$ = newStruct($1->name);
          //    $$->level = NestLevel;
          //    $$->tagsym = $1;
          //    //addSym (StructTab, $$, $$->tag, $$->level, currBlockno, 0);
          //  }
        }
#line 6158 "y.tab.c"
    break;

  case 376: /* opt_assign_expr: '=' constant_expr  */
#line 2498 "SDCC.y"
        {
          value *val;

          val = constExprValue((yyvsp[0].asts), true);
          if (!val) // Not a constant expression
            {
              werror (E_CONST_EXPECTED);
              val = constIntVal("0");
            }
          else if (!IS_INT(val->type) && !IS_CHAR(val->type) && !IS_BOOL(val->type))
            {
              werror(E_ENUM_NON_INTEGER);
              SNPRINTF(lbuff, sizeof(lbuff), "%lld", (long long int) ullFromVal(val));
              val = constVal(lbuff);
            }
          (yyval.val) = cenum = val;
        }
#line 6180 "y.tab.c"
    break;

  case 377: /* opt_assign_expr: %empty  */
#line 2515 "SDCC.y"
        {
          if (cenum)
            {
              SNPRINTF(lbuff, sizeof(lbuff), "%lld", (long long int) ullFromVal(cenum)+1);
              (yyval.val) = cenum = constVal(lbuff);
            }
          else
            {
              (yyval.val) = cenum = constCharVal(0);
            }
        }
#line 6196 "y.tab.c"
    break;

  case 378: /* specifier_qualifier_list: type_specifier_list_  */
#line 2528 "SDCC.y"
                                                { (yyval.lnk) = finalizeSpec((yyvsp[0].lnk)); }
#line 6202 "y.tab.c"
    break;

  case 380: /* type_specifier_list_: type_specifier_list_ type_specifier_qualifier  */
#line 2533 "SDCC.y"
                                                   {
     /* if the decl $2 is not a specifier */
     /* find the spec and replace it      */
     (yyval.lnk) = mergeDeclSpec((yyvsp[-1].lnk), (yyvsp[0].lnk), "type_specifier_list type_specifier skipped");
   }
#line 6212 "y.tab.c"
    break;

  case 382: /* identifier_list: identifier_list ',' identifier  */
#line 2543 "SDCC.y"
         {
           (yyvsp[0].sym)->next = (yyvsp[-2].sym);
           (yyval.sym) = (yyvsp[0].sym);
         }
#line 6221 "y.tab.c"
    break;

  case 383: /* type_name: declaration_specifiers  */
#line 2551 "SDCC.y"
        {
          if (IS_SPEC ((yyvsp[0].lnk)) && !IS_VALID_PARAMETER_STORAGE_CLASS_SPEC ((yyvsp[0].lnk)))
            {
              werror (E_STORAGE_CLASS_FOR_PARAMETER, "type name");
            }
          (yyval.lnk) = (yyvsp[0].lnk); ignoreTypedefType = 0;
        }
#line 6233 "y.tab.c"
    break;

  case 384: /* type_name: declaration_specifiers abstract_declarator  */
#line 2559 "SDCC.y"
        {
          /* go to the end of the list */
          sym_link *p;

          if (IS_SPEC ((yyvsp[-1].lnk)) && !IS_VALID_PARAMETER_STORAGE_CLASS_SPEC ((yyvsp[-1].lnk)))
            {
              werror (E_STORAGE_CLASS_FOR_PARAMETER, "type name");
            }
          pointerTypes ((yyvsp[0].lnk),(yyvsp[-1].lnk));
          for (p = (yyvsp[0].lnk); p && p->next; p = p->next)
            ;
          if (!p)
            {
              werror(E_SYNTAX_ERROR, yytext);
            }
          else
            {
              p->next = (yyvsp[-1].lnk);
            }
          (yyval.lnk) = (yyvsp[0].lnk);
          ignoreTypedefType = 0;
        }
#line 6260 "y.tab.c"
    break;

  case 385: /* critical: CRITICAL  */
#line 2584 "SDCC.y"
                {
                   if (inCriticalFunction || inCriticalBlock)
                     werror(E_INVALID_CRITICAL);
                   inCriticalBlock = 1;
                   STACK_PUSH(continueStack,NULL);
                   STACK_PUSH(breakStack,NULL);
                   (yyval.sym) = NULL;
                }
#line 6273 "y.tab.c"
    break;

  case 386: /* critical_statement: critical statement  */
#line 2595 "SDCC.y"
                         {
                   STACK_POP(breakStack);
                   STACK_POP(continueStack);
                   (yyval.asts) = newNode(CRITICAL,(yyvsp[0].asts),NULL);
                   inCriticalBlock = 0;
                }
#line 6284 "y.tab.c"
    break;

  case 388: /* statements_and_implicit: statement_list implicit_block  */
#line 2606 "SDCC.y"
     {
       (yyval.asts) = newNode(NULLOP, (yyvsp[-1].asts), (yyvsp[0].asts));
       if (!options.std_c99)
         werror(E_DECL_AFTER_STATEMENT_C99);
     }
#line 6294 "y.tab.c"
    break;

  case 389: /* $@18: %empty  */
#line 2614 "SDCC.y"
     {
       NestLevel += SUBLEVEL_UNIT;
       STACK_PUSH(blockNum, currBlockno);
       btree_add_child(currBlockno, ++blockNo);
       currBlockno = blockNo;
       ignoreTypedefType = 0;
     }
#line 6306 "y.tab.c"
    break;

  case 390: /* declaration_after_statement: $@18 declaration_list  */
#line 2621 "SDCC.y"
                                              { (yyval.sym) = (yyvsp[0].sym); }
#line 6312 "y.tab.c"
    break;

  case 391: /* implicit_block: declaration_after_statement statements_and_implicit  */
#line 2626 "SDCC.y"
     {
       leaveBlockScope (currBlockno);
       NestLevel -= SUBLEVEL_UNIT;
       currBlockno = STACK_POP(blockNum);
       (yyval.asts) = createBlock((yyvsp[-1].sym), (yyvsp[0].asts));
       cleanUpLevel(StructTab, NestLevel + SUBLEVEL_UNIT);
     }
#line 6324 "y.tab.c"
    break;

  case 392: /* implicit_block: declaration_after_statement  */
#line 2634 "SDCC.y"
     {
       leaveBlockScope (currBlockno);
       NestLevel -= SUBLEVEL_UNIT;
       currBlockno = STACK_POP(blockNum);
       (yyval.asts) = createBlock((yyvsp[0].sym), NULL);
       cleanUpLevel(StructTab, NestLevel + SUBLEVEL_UNIT);
     }
#line 6336 "y.tab.c"
    break;

  case 393: /* declaration_list: declaration  */
#line 2645 "SDCC.y"
     {
       /* if this is typedef declare it immediately */
       if ( (yyvsp[0].sym) && IS_TYPEDEF((yyvsp[0].sym)->etype)) {
         allocVariables ((yyvsp[0].sym));
         (yyval.sym) = NULL;
       }
       else
         (yyval.sym) = (yyvsp[0].sym);
       ignoreTypedefType = 0;/*printf("1 %s %d %d\n", $1->name, IS_ARRAY ($1->type), DCL_ARRAY_VLA ($1->type));
       if (IS_ARRAY ($1->type) && DCL_ARRAY_VLA ($1->type))
         werror (E_VLA_OBJECT);
       else*/
         addSymChain(&(yyvsp[0].sym));
     }
#line 6355 "y.tab.c"
    break;

  case 394: /* declaration_list: declaration_list declaration  */
#line 2661 "SDCC.y"
     {
       symbol   *sym;

       /* if this is a typedef */
       if ((yyvsp[0].sym) && IS_TYPEDEF((yyvsp[0].sym)->etype)) {
         allocVariables ((yyvsp[0].sym));
         (yyval.sym) = (yyvsp[-1].sym);
       }
       else {
         /* get to the end of the previous decl */
         if ( (yyvsp[-1].sym) ) {
           (yyval.sym) = sym = (yyvsp[-1].sym);
           while (sym->next)
             sym = sym->next;
           sym->next = (yyvsp[0].sym);
         }
         else
           (yyval.sym) = (yyvsp[0].sym);
       }
       ignoreTypedefType = 0;/*printf("2 %s %d %d\n", $1->name, IS_ARRAY ($1->type), DCL_ARRAY_VLA ($1->type));
       if (IS_ARRAY ($2->type) && DCL_ARRAY_VLA ($2->type))
         werror (E_VLA_OBJECT);
       else*/
         addSymChain(&(yyvsp[0].sym));
     }
#line 6385 "y.tab.c"
    break;

  case 395: /* kr_declaration: declaration_specifiers init_declarator_list ';'  */
#line 2695 "SDCC.y"
      {
         /* add the specifier list to the id */
         symbol *sym , *sym1;

         for (sym1 = sym = reverseSyms((yyvsp[-1].sym));sym != NULL;sym = sym->next) {
             sym_link *lnk = copyLinkChain((yyvsp[-2].lnk));
             sym_link *l0 = NULL, *l1 = NULL, *l2 = NULL;
             /* check illegal declaration */
             for (l0 = sym->type; l0 != NULL; l0 = l0->next)
               if (IS_PTR (l0))
                 break;
             /* check if creating instances of structs with flexible arrays */
             for (l1 = lnk; l1 != NULL; l1 = l1->next)
               if (IS_STRUCT (l1) && SPEC_STRUCT (l1)->b_flexArrayMember)
                 break;
             if (!options.std_c99 && l0 == NULL && l1 != NULL && SPEC_EXTR((yyvsp[-2].lnk)) != 1)
               werror (W_FLEXARRAY_INSTRUCT, sym->name);
             /* check if creating instances of function type */
             for (l1 = lnk; l1 != NULL; l1 = l1->next)
               if (IS_FUNC (l1))
                 break;
             for (l2 = lnk; l2 != NULL; l2 = l2->next)
               if (IS_PTR (l2))
                 break;
             if (l0 == NULL && l2 == NULL && l1 != NULL)
               werrorfl(sym->fileDef, sym->lineDef, E_TYPE_IS_FUNCTION, sym->name);
             /* do the pointer stuff */
             pointerTypes(sym->type,lnk);
             addDecl (sym,0,lnk);
         }

         uselessDecl = true;
         (yyval.sym) = sym1;
      }
#line 6424 "y.tab.c"
    break;

  case 396: /* kr_declaration_list: kr_declaration  */
#line 2733 "SDCC.y"
     {
       /* if this is typedef declare it immediately */
       if ( (yyvsp[0].sym) && IS_TYPEDEF((yyvsp[0].sym)->etype)) {
         allocVariables ((yyvsp[0].sym));
         (yyval.sym) = NULL;
       }
       else
         (yyval.sym) = (yyvsp[0].sym);
       ignoreTypedefType = 0;
       addSymChain(&(yyvsp[0].sym));
     }
#line 6440 "y.tab.c"
    break;

  case 397: /* kr_declaration_list: kr_declaration_list kr_declaration  */
#line 2745 "SDCC.y"
     {
       symbol   *sym;

       /* if this is a typedef */
       if ((yyvsp[0].sym) && IS_TYPEDEF((yyvsp[0].sym)->etype)) {
         allocVariables ((yyvsp[0].sym));
         (yyval.sym) = (yyvsp[-1].sym);
       }
       else {
         /* get to the end of the previous decl */
         if ( (yyvsp[-1].sym) ) {
           (yyval.sym) = sym = (yyvsp[-1].sym);
           while (sym->next)
             sym = sym->next;
           sym->next = (yyvsp[0].sym);
         }
         else
           (yyval.sym) = (yyvsp[0].sym);
       }
       ignoreTypedefType = 0;
       addSymChain(&(yyvsp[0].sym));
     }
#line 6467 "y.tab.c"
    break;

  case 398: /* statement_list: unlabeled_statement  */
#line 2770 "SDCC.y"
                                        { (yyval.asts) = (yyvsp[0].asts); }
#line 6473 "y.tab.c"
    break;

  case 399: /* statement_list: label  */
#line 2771 "SDCC.y"
                                        { (yyval.asts) = (yyvsp[0].asts); }
#line 6479 "y.tab.c"
    break;

  case 400: /* statement_list: statement_list unlabeled_statement  */
#line 2772 "SDCC.y"
                                        { (yyval.asts) = newNode(NULLOP,(yyvsp[-1].asts),(yyvsp[0].asts));}
#line 6485 "y.tab.c"
    break;

  case 401: /* statement_list: statement_list label  */
#line 2773 "SDCC.y"
                                        { (yyval.asts) = newNode(NULLOP,(yyvsp[-1].asts),(yyvsp[0].asts));}
#line 6491 "y.tab.c"
    break;

  case 402: /* while: WHILE  */
#line 2776 "SDCC.y"
               {  /* create and push the continue , break & body labels */
                  static int Lblnum = 0;
                  /* continue */
                  SNPRINTF (lbuff, sizeof(lbuff), "_whilecontinue_%d",Lblnum);
                  STACK_PUSH(continueStack,newSymbol(lbuff,NestLevel));
                  /* break */
                  SNPRINTF (lbuff, sizeof(lbuff), "_whilebreak_%d",Lblnum);
                  STACK_PUSH(breakStack,newSymbol(lbuff,NestLevel));
                  /* body */
                  SNPRINTF (lbuff, sizeof(lbuff), "_whilebody_%d",Lblnum++);
                  (yyval.sym) = newSymbol(lbuff,NestLevel);
               }
#line 6508 "y.tab.c"
    break;

  case 403: /* do: DO  */
#line 2790 "SDCC.y"
        {  /* create and push the continue , break & body Labels */
           static int Lblnum = 0;

           /* continue */
           SNPRINTF(lbuff, sizeof(lbuff), "_docontinue_%d",Lblnum);
           STACK_PUSH(continueStack,newSymbol(lbuff,NestLevel));
           /* break */
           SNPRINTF(lbuff, sizeof(lbuff), "_dobreak_%d",Lblnum);
           STACK_PUSH(breakStack,newSymbol(lbuff,NestLevel));
           /* do body */
           SNPRINTF(lbuff, sizeof(lbuff), "_dobody_%d",Lblnum++);
           (yyval.sym) = newSymbol (lbuff,NestLevel);
        }
#line 6526 "y.tab.c"
    break;

  case 404: /* for: FOR  */
#line 2805 "SDCC.y"
          { /* create & push continue, break & body labels */
            static int Lblnum = 0;

           NestLevel += LEVEL_UNIT;
           STACK_PUSH(blockNum, currBlockno);
           btree_add_child(currBlockno, ++blockNo);
           currBlockno = blockNo;
           ignoreTypedefType = 0;

            /* continue */
            SNPRINTF(lbuff, sizeof(lbuff), "_forcontinue_%d",Lblnum);
            STACK_PUSH(continueStack,newSymbol(lbuff,NestLevel));
            /* break    */
            SNPRINTF(lbuff, sizeof(lbuff), "_forbreak_%d",Lblnum);
            STACK_PUSH(breakStack,newSymbol(lbuff,NestLevel));
            /* body */
            SNPRINTF(lbuff, sizeof(lbuff), "_forbody_%d",Lblnum);
            (yyval.sym) = newSymbol(lbuff,NestLevel);
            /* condition */
            SNPRINTF(lbuff, sizeof(lbuff), "_forcond_%d",Lblnum++);
            STACK_PUSH(forStack,newSymbol(lbuff,NestLevel));
          }
#line 6553 "y.tab.c"
    break;

  case 406: /* asm_statement: ASM '(' asm_string_literal ')' ';'  */
#line 2835 "SDCC.y"
      {
        ast *ex;

        seqPointNo++;
        ex = newNode(INLINEASM, NULL, NULL);
        ex->values.inlineasm = strdup(copyStr ((yyvsp[-2].yystr), NULL));
        seqPointNo++;
        (yyval.asts) = ex;
     }
#line 6567 "y.tab.c"
    break;

  case 407: /* asm_statement: INLINEASM ';'  */
#line 2845 "SDCC.y"
      {
        ast *ex;

        seqPointNo++;
        ex = newNode(INLINEASM, NULL, NULL);
        ex->values.inlineasm = strdup((yyvsp[-1].yystr));
        seqPointNo++;
        (yyval.asts) = ex;
      }
#line 6581 "y.tab.c"
    break;

  case 408: /* addressmod: ADDRESSMOD identifier identifier ';'  */
#line 2857 "SDCC.y"
                                          {
     symbol *sym;
     if ((sym = findSymWithLevel (AddrspaceTab, (yyvsp[-1].sym))) && sym->level == (yyvsp[-1].sym)->level)
       werrorfl (sym->fileDef, sym->lineDef, E_PREVIOUS_DEF);
     if (!findSymWithLevel (SymbolTab, (yyvsp[-2].sym)))
       werror (E_ID_UNDEF, (yyvsp[-2].sym)->name);
     addSym (AddrspaceTab, (yyvsp[-1].sym), (yyvsp[-1].sym)->name, (yyvsp[-1].sym)->level, (yyvsp[-1].sym)->block, 0);
     sym = findSymWithLevel (AddrspaceTab, (yyvsp[-1].sym));
     sym->addressmod[0] = findSymWithLevel (SymbolTab, (yyvsp[-2].sym));
     sym->addressmod[1] = 0;
   }
#line 6597 "y.tab.c"
    break;

  case 409: /* addressmod: ADDRESSMOD identifier SD_CONST identifier ';'  */
#line 2868 "SDCC.y"
                                                   {
     symbol *sym;
     sym_link *type;
     if ((sym = findSymWithLevel (AddrspaceTab, (yyvsp[-1].sym))) && sym->level == (yyvsp[-1].sym)->level)
       werrorfl (sym->fileDef, sym->lineDef, E_PREVIOUS_DEF);
     if (!findSymWithLevel (SymbolTab, (yyvsp[-3].sym)))
       werror (E_ID_UNDEF, (yyvsp[-3].sym)->name);
     addSym (AddrspaceTab, (yyvsp[-1].sym), (yyvsp[-1].sym)->name, (yyvsp[-1].sym)->level, (yyvsp[-1].sym)->block, 0);
     sym = findSymWithLevel (AddrspaceTab, (yyvsp[-1].sym));
     sym->addressmod[0] = findSymWithLevel (SymbolTab, (yyvsp[-3].sym));
     sym->addressmod[1] = 0;
     type = newLink (SPECIFIER);
     SPEC_CONST(type) = 1;
     sym->type = sym->etype = type;
   }
#line 6617 "y.tab.c"
    break;

  case 410: /* identifier: IDENTIFIER  */
#line 2886 "SDCC.y"
                  { (yyval.sym) = newSymbol ((yyvsp[0].yychar), NestLevel); }
#line 6623 "y.tab.c"
    break;


#line 6627 "y.tab.c"

      default: break;
    }
  /* User semantic actions sometimes alter yychar, and that requires
     that yytoken be updated with the new translation.  We take the
     approach of translating immediately before every use of yytoken.
     One alternative is translating here after every semantic action,
     but that translation would be missed if the semantic action invokes
     YYABORT, YYACCEPT, or YYERROR immediately after altering yychar or
     if it invokes YYBACKUP.  In the case of YYABORT or YYACCEPT, an
     incorrect destructor might then be invoked immediately.  In the
     case of YYERROR or YYBACKUP, subsequent parser actions might lead
     to an incorrect destructor call or verbose syntax error message
     before the lookahead is translated.  */
  YY_SYMBOL_PRINT ("-> $$ =", YY_CAST (yysymbol_kind_t, yyr1[yyn]), &yyval, &yyloc);

  YYPOPSTACK (yylen);
  yylen = 0;

  *++yyvsp = yyval;

  /* Now 'shift' the result of the reduction.  Determine what state
     that goes to, based on the state we popped back to and the rule
     number reduced by.  */
  {
    const int yylhs = yyr1[yyn] - YYNTOKENS;
    const int yyi = yypgoto[yylhs] + *yyssp;
    yystate = (0 <= yyi && yyi <= YYLAST && yycheck[yyi] == *yyssp
               ? yytable[yyi]
               : yydefgoto[yylhs]);
  }

  goto yynewstate;


/*--------------------------------------.
| yyerrlab -- here on detecting error.  |
`--------------------------------------*/
yyerrlab:
  /* Make sure we have latest lookahead translation.  See comments at
     user semantic actions for why this is necessary.  */
  yytoken = yychar == YYEMPTY ? YYSYMBOL_YYEMPTY : YYTRANSLATE (yychar);
  /* If not already recovering from an error, report this error.  */
  if (!yyerrstatus)
    {
      ++yynerrs;
      yyerror (YY_("syntax error"));
    }

  if (yyerrstatus == 3)
    {
      /* If just tried and failed to reuse lookahead token after an
         error, discard it.  */

      if (yychar <= YYEOF)
        {
          /* Return failure if at end of input.  */
          if (yychar == YYEOF)
            YYABORT;
        }
      else
        {
          yydestruct ("Error: discarding",
                      yytoken, &yylval);
          yychar = YYEMPTY;
        }
    }

  /* Else will try to reuse lookahead token after shifting the error
     token.  */
  goto yyerrlab1;


/*---------------------------------------------------.
| yyerrorlab -- error raised explicitly by YYERROR.  |
`---------------------------------------------------*/
yyerrorlab:
  /* Pacify compilers when the user code never invokes YYERROR and the
     label yyerrorlab therefore never appears in user code.  */
  if (0)
    YYERROR;
  ++yynerrs;

  /* Do not reclaim the symbols of the rule whose action triggered
     this YYERROR.  */
  YYPOPSTACK (yylen);
  yylen = 0;
  YY_STACK_PRINT (yyss, yyssp);
  yystate = *yyssp;
  goto yyerrlab1;


/*-------------------------------------------------------------.
| yyerrlab1 -- common code for both syntax error and YYERROR.  |
`-------------------------------------------------------------*/
yyerrlab1:
  yyerrstatus = 3;      /* Each real token shifted decrements this.  */

  /* Pop stack until we find a state that shifts the error token.  */
  for (;;)
    {
      yyn = yypact[yystate];
      if (!yypact_value_is_default (yyn))
        {
          yyn += YYSYMBOL_YYerror;
          if (0 <= yyn && yyn <= YYLAST && yycheck[yyn] == YYSYMBOL_YYerror)
            {
              yyn = yytable[yyn];
              if (0 < yyn)
                break;
            }
        }

      /* Pop the current state because it cannot handle the error token.  */
      if (yyssp == yyss)
        YYABORT;


      yydestruct ("Error: popping",
                  YY_ACCESSING_SYMBOL (yystate), yyvsp);
      YYPOPSTACK (1);
      yystate = *yyssp;
      YY_STACK_PRINT (yyss, yyssp);
    }

  YY_IGNORE_MAYBE_UNINITIALIZED_BEGIN
  *++yyvsp = yylval;
  YY_IGNORE_MAYBE_UNINITIALIZED_END


  /* Shift the error token.  */
  YY_SYMBOL_PRINT ("Shifting", YY_ACCESSING_SYMBOL (yyn), yyvsp, yylsp);

  yystate = yyn;
  goto yynewstate;


/*-------------------------------------.
| yyacceptlab -- YYACCEPT comes here.  |
`-------------------------------------*/
yyacceptlab:
  yyresult = 0;
  goto yyreturnlab;


/*-----------------------------------.
| yyabortlab -- YYABORT comes here.  |
`-----------------------------------*/
yyabortlab:
  yyresult = 1;
  goto yyreturnlab;


/*-----------------------------------------------------------.
| yyexhaustedlab -- YYNOMEM (memory exhaustion) comes here.  |
`-----------------------------------------------------------*/
yyexhaustedlab:
  yyerror (YY_("memory exhausted"));
  yyresult = 2;
  goto yyreturnlab;


/*----------------------------------------------------------.
| yyreturnlab -- parsing is finished, clean up and return.  |
`----------------------------------------------------------*/
yyreturnlab:
  if (yychar != YYEMPTY)
    {
      /* Make sure we have latest lookahead translation.  See comments at
         user semantic actions for why this is necessary.  */
      yytoken = YYTRANSLATE (yychar);
      yydestruct ("Cleanup: discarding lookahead",
                  yytoken, &yylval);
    }
  /* Do not reclaim the symbols of the rule whose action triggered
     this YYABORT or YYACCEPT.  */
  YYPOPSTACK (yylen);
  YY_STACK_PRINT (yyss, yyssp);
  while (yyssp != yyss)
    {
      yydestruct ("Cleanup: popping",
                  YY_ACCESSING_SYMBOL (+*yyssp), yyvsp);
      YYPOPSTACK (1);
    }
#ifndef yyoverflow
  if (yyss != yyssa)
    YYSTACK_FREE (yyss);
#endif

  return yyresult;
}

#line 2888 "SDCC.y"

