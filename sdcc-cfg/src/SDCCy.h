/* A Bison parser, made by GNU Bison 3.8.2.  */

/* Bison interface for Yacc-like parsers in C

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

/* DO NOT RELY ON FEATURES THAT ARE NOT DOCUMENTED in the manual,
   especially those whose name start with YY_ or yy_.  They are
   private implementation details that can be changed or removed.  */

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

#line 401 "y.tab.h"

};
typedef union YYSTYPE YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define YYSTYPE_IS_DECLARED 1
#endif


extern YYSTYPE yylval;


int yyparse (void);


#endif /* !YY_YY_Y_TAB_H_INCLUDED  */
