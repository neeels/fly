
const uint8_t font[][5] = {
  { 0x00,0x00,0x00,0x00,0x00 }, //0/ -->
  { 0x00,0x00,0x00,0x00,0x00 }, //1/ --> ^A
  { 0x00,0x00,0x00,0x00,0x00 }, //2/ --> ^B
  { 0x00,0x00,0x00,0x00,0x00 }, //3/ --> ^C
  { 0x00,0x00,0x00,0x00,0x00 }, //4/ --> ^D
  { 0x00,0x00,0x00,0x00,0x00 }, //5/ --> ^E
  { 0x00,0x00,0x00,0x00,0x00 }, //6/ --> ^F
  { 0x00,0x00,0x00,0x00,0x00 }, //7/ --> ^G
  { 0x00,0x00,0x00,0x00,0x00 }, //8/ --> ^H
  { 0x00,0x00,0x00,0x00,0x00 }, //9/ -->
  { 0x00,0x00,0x00,0x00,0x00 }, //10/ -->
  { 0x00,0x00,0x00,0x00,0x00 }, //11/ -->
  { 0x00,0x00,0x00,0x00,0x00 }, //12/ -->
  { 0x00,0x00,0x00,0x00,0x00 }, //13/ -->  
  { 0x00,0x00,0x00,0x00,0x00 }, //14/ --> ^N
  { 0x00,0x00,0x00,0x00,0x00 }, //15/ --> ^O
  { 0x3f,0x3f,0x30,0x3f,0x3f }, //16/ --> u
  { 0xff,0xff,0x33,0x3f,0x3f }, //17/ --> p
  { 0xfc,0xfc,0xcc,0xff,0xff }, //18/ --> d
  { 0xfc,0xfc,0x0c,0xfc,0xfc }, //19/ --> n
  { 0x20,0x20,0x20,0x24,0x2A }, //20/ --> ^T
  { 0x20,0x20,0x24,0x2A,0x11 }, //21/ --> ^U
  { 0x20,0x20,0x20,0x24,0x2A }, //22/ --> ^V
  { 0x20,0x30,0x20,0x30,0x28 }, //23/ --> ^W
  { 0x24,0x22,0x21,0x24,0x2A }, //24/ --> ^X
  { 0x22,0x21,0x24,0x2A,0x11 }, //25/ --> ^Y
  { 0x24,0x22,0x21,0x24,0x2A }, //26/ --> ^Z
  { 0x80,0x80,0x40,0x30,0x20 }, //27/ --> ^[
  { 0x20,0x20,0x60,0xA0,0x60 }, //28/ -->  
  { 0x20,0x20,0x60,0xB0,0x60 }, //29/ --> ^]
  { 0x30,0x28,0x60,0xA0,0x60 }, //30/ --> ^^
  { 0x04,0x06,0x1D,0x25,0x24 }, //31/ --> ^_
  { 0x00,0x00,0x00,0x00,0x00 }, //32/ -->
  { 0x00,0x4F,0x00,0x00,0x00 }, //33/ --> !
  { 0x07,0x00,0x07,0x00,0x00 }, //34/ --> "
  { 0x14,0x7F,0x14,0x7F,0x14 }, //35/ --> #
  { 0x24,0x2A,0x7F,0x2A,0x12 }, //36/ --> $
  { 0x22,0x10,0x08,0x04,0x22 }, //37/ --> %
  { 0x36,0x49,0x55,0x22,0x40 }, //38/ --> &
  { 0x00,0x05,0x03,0x00,0x00 }, //39/ --> '
  { 0x1C,0x22,0x41,0x00,0x00 }, //40/ --> (
  { 0x41,0x22,0x1C,0x00,0x00 }, //41/ --> )
  { 0x14,0x08,0x3E,0x08,0x14 }, //42/ --> *
  { 0x08,0x08,0x3E,0x08,0x08 }, //43/ --> +
  { 0xa0,0x60,0x00,0x00,0x00 }, //44/ --> ,
  { 0x08,0x08,0x08,0x08,0x08 }, //45/ --> -
  { 0x80,0x00,0x00,0x00,0x00 }, //46/ --> .
  { 0x20,0x10,0x08,0x04,0x02 }, //47/ --> /
  { 0x3E,0x51,0x49,0x45,0x3E }, //48/ --> 0
  { 0x00,0x42,0x7F,0x40,0x00 }, //49/ --> 1
  { 0x42,0x61,0x51,0x49,0x46 }, //50/ --> 2
  { 0x21,0x41,0x45,0x4B,0x31 }, //51/ --> 3
  { 0x18,0x14,0x12,0x7F,0x10 }, //52/ --> 4
  { 0x27,0x45,0x45,0x45,0x39 }, //53/ --> 5
  { 0x3C,0x4A,0x49,0x49,0x30 }, //54/ --> 6
  { 0x01,0x71,0x09,0x05,0x03 }, //55/ --> 7
  { 0x36,0x49,0x49,0x49,0x36 }, //56/ --> 8
  { 0x06,0x49,0x49,0x29,0x1E }, //57/ --> 9
  { 0x00,0x36,0x36,0x00,0x00 }, //58/ --> :
  { 0x00,0x56,0x36,0x00,0x00 }, //59/ --> ;
  { 0x08,0x14,0x22,0x41,0x00 }, //60/ --> <
  { 0x24,0x24,0x24,0x24,0x24 }, //61/ --> =
  { 0x00,0x41,0x22,0x14,0x08 }, //62/ --> >
  { 0x02,0x01,0x51,0x09,0x06 }, //63/ --> ?
  { 0x32,0x49,0x79,0x41,0x3E }, //64/ --> @
  { 0x7E,0x11,0x11,0x11,0x7E }, //65/ --> A
  { 0x7F,0x49,0x49,0x49,0x36 }, //66/ --> B
  { 0x3E,0x41,0x41,0x41,0x22 }, //67/ --> C
  { 0x7F,0x41,0x41,0x22,0x1C }, //68/ --> D
  { 0x7F,0x49,0x49,0x49,0x41 }, //69/ --> E
  { 0x7F,0x09,0x09,0x09,0x01 }, //70/ --> F
  { 0x3E,0x41,0x49,0x49,0x3A }, //71/ --> G
  { 0x7F,0x08,0x08,0x08,0x7F }, //72/ --> H
  { 0x00,0x41,0x7F,0x41,0x00 }, //73/ --> I
  { 0x20,0x40,0x41,0x3F,0x01 }, //74/ --> J
  { 0x7F,0x08,0x14,0x22,0x41 }, //75/ --> K
  { 0x7F,0x40,0x40,0x40,0x40 }, //76/ --> L
  { 0x7F,0x02,0x0C,0x02,0x7F }, //77/ --> M
  { 0x7F,0x04,0x08,0x10,0x7F }, //78/ --> N
  { 0x3E,0x41,0x41,0x41,0x3E }, //79/ --> O
  { 0x7F,0x09,0x09,0x09,0x06 }, //80/ --> P
  { 0x3E,0x41,0x51,0x21,0x5E }, //81/ --> Q
  { 0x7F,0x09,0x19,0x29,0x46 }, //82/ --> R
  { 0x46,0x49,0x49,0x49,0x31 }, //83/ --> S
  { 0x01,0x01,0x7F,0x01,0x01 }, //84/ --> T
  { 0x3F,0x40,0x40,0x40,0x3F }, //85/ --> U
  { 0x1F,0x20,0x40,0x20,0x1F }, //86/ --> V
  { 0x3F,0x40,0x60,0x40,0x3F }, //87/ --> W
  { 0x63,0x14,0x08,0x14,0x63 }, //88/ --> X
  { 0x07,0x08,0x70,0x08,0x07 }, //89/ --> Y
  { 0x61,0x51,0x49,0x45,0x43 }, //90/ --> Z
  { 0x7F,0x41,0x41,0x00,0x00 }, //91/ --> [
  { 0x02,0x04,0x08,0x10,0x20 }, //92/ --> '\'
  { 0x41,0x41,0x7F,0x00,0x00 }, //93/ --> ]
  { 0x04,0x02,0x01,0x02,0x04 }, //94/ --> ^
  { 0x40,0x40,0x40,0x40,0x40 }, //95/ --> _
  { 0x01,0x02,0x04,0x00,0x00 }, //96/ --> `
  { 0x20,0x54,0x54,0x54,0x78 }, //97/ --> a
  { 0x7F,0x44,0x44,0x44,0x38 }, //98/ --> b
  { 0x38,0x44,0x44,0x44,0x00 }, //99/ --> c
  { 0x38,0x44,0x44,0x48,0x7F }, //100/ --> d
  { 0x38,0x54,0x54,0x54,0x18 }, //101/ --> e
  { 0x10,0x7E,0x11,0x01,0x02 }, //102/ --> f
  { 0x0C,0x52,0x52,0x52,0x3E }, //103/ --> g
  { 0x7F,0x08,0x04,0x04,0x78 }, //104/ --> h
  { 0x00,0x44,0x7D,0x40,0x00 }, //105/ --> i
  { 0x20,0x40,0x40,0x3D,0x00 }, //106/ --> j
  { 0x7F,0x10,0x28,0x44,0x00 }, //107/ --> k
  { 0x00,0x41,0x7F,0x40,0x00 }, //108/ --> l
  { 0x7C,0x04,0x18,0x04,0x78 }, //109/ --> m
  { 0x7C,0x08,0x04,0x04,0x78 }, //110/ --> n
  { 0x38,0x44,0x44,0x44,0x38 }, //111/ --> o
  { 0x7C,0x14,0x14,0x14,0x08 }, //112/ --> p
  { 0x08,0x14,0x14,0x18,0x7C }, //113/ --> q
  { 0x7C,0x08,0x04,0x04,0x08 }, //114/ --> r
  { 0x48,0x54,0x54,0x54,0x20 }, //115/ --> s
  { 0x04,0x3F,0x44,0x40,0x20 }, //116/ --> t
  { 0x3C,0x40,0x40,0x20,0x7C }, //117/ --> u
  { 0x1C,0x20,0x40,0x20,0x1C }, //118/ --> v
  { 0x3c,0x40,0x20,0x40,0x3c }, //119/ --> w
  { 0x44,0x28,0x10,0x28,0x44 }, //120/ --> x
  { 0x06,0x48,0x48,0x48,0x3E }, //121/ --> y
  { 0x44,0x64,0x54,0x4C,0x44 }, //122/ --> z
  { 0x08,0x36,0x41,0x00,0x00 }, //123/ --> {
  { 0x00,0x7F,0x00,0x00,0x00 }, //124/ --> |
  { 0x41,0x36,0x08,0x00,0x00 }, //125/ --> }
  { 0x08,0x08,0x2A,0x1C,0x08 }, //126/ --> ~
  { 0x08,0x1C,0x2A,0x08,0x08 }, //127/ --> ^?
/*
  { 0x3C,0x42,0x41,0x42,0x3C }, //128/ --> <80>
  { 0x30,0x28,0x60,0xA0,0x60 }, //129/ --> <81>
  { 0x20,0x20,0x20,0xA0,0x20 }, //130/ --> <82>
  { 0x20,0x20,0x20,0xB0,0x20 }, //131/ --> <83>
  { 0x30,0x28,0x20,0xA0,0x20 }, //132/ --> <84>
  { 0x20,0x20,0x22,0x20,0x22 }, //133/ --> <85>
  { 0x20,0x20,0x22,0x30,0x22 }, //134/ --> <86>
  { 0x30,0x28,0x22,0x20,0x22 }, //135/ --> <87>
  { 0x20,0x20,0x22,0x21,0x22 }, //136/ --> <88>
  { 0x20,0x20,0x22,0x31,0x22 }, //137/ --> <89>
  { 0x30,0x28,0x22,0x21,0x22 }, //138/ --> <8a>
  { 0x20,0x28,0x28,0x28,0xB0 }, //139/ --> <8b>
  { 0x20,0x28,0x28,0x28,0x30 }, //140/ --> <8c>
  { 0xC0,0xA8,0x28,0x68,0xB0 }, //141/ --> <8d>
  { 0x00,0x80,0x80,0x44,0x32 }, //142/ --> <8e>
  { 0x24,0x25,0x24,0x38,0x20 }, //143/ --> <8f>
  { 0x24,0x22,0x21,0x24,0x2A }, //144/ --> <90>
  { 0x80,0x80,0x40,0x34,0x20 }, //145/ --> <91>
  { 0x20,0x20,0x38,0x20,0x38 }, //146/ --> <92>
  { 0x20,0x38,0x20,0x38,0x20 }, //147/ --> <93>
  { 0x80,0x80,0x78,0x20,0x38 }, //148/ --> <94>
  { 0x20,0x20,0x38,0x22,0x39 }, //149/ --> <95>
  { 0x20,0x38,0x22,0x39,0x22 }, //150/ --> <96>
  { 0x80,0x80,0x78,0x22,0x39 }, //151/ --> <97>
  { 0x20,0x20,0x20,0x24,0x2A }, //152/ --> <98>
  { 0x30,0x20,0x30,0x28,0x28 }, //153/ --> <99>
  { 0x80,0x80,0x60,0x30,0x28 }, //154/ --> <9a>
  { 0x20,0x30,0x20,0x30,0x28 }, //155/ --> <9b>
  { 0x30,0x20,0x30,0x28,0x2A }, //156/ --> <9c>
  { 0x80,0x80,0x60,0x30,0x28 }, //157/ --> <9d>
  { 0x20,0x3E,0x30,0x28,0x28 }, //158/ --> <9e>
  { 0x20,0x3E,0x30,0x28,0x2A }, //159/ --> <9f>
  { 0x20,0x20,0x20,0x30,0x28 }, //160/ -->
  { 0x20,0x20,0x30,0x28,0x28 }, //161/ --> ¡
  { 0x40,0xA0,0xB0,0x28,0x28 }, //162/ --> ¢
  { 0x20,0x20,0x20,0x30,0x28 }, //163/ --> £
  { 0x20,0x20,0x30,0x28,0x2A }, //164/ --> ¤
  { 0x40,0xA0,0xB0,0x28,0x2A }, //165/ --> ¥
  { 0x20,0x20,0x20,0x30,0x28 }, //166/ --> ¦
  { 0x20,0x30,0x28,0x2A,0x30 }, //167/ --> §
  { 0x18,0x20,0x20,0x30,0x28 }, //168/ --> ¨
  { 0x20,0x20,0x20,0x30,0x2A }, //169/ --> ©
  { 0x20,0x30,0x2A,0x28,0x32 }, //170/ --> ª
  { 0x80,0x80,0xB2,0xA8,0x7A }, //171/ --> «
  { 0x25,0x25,0x25,0x25,0x25 }, //172/ --> ¬
  { 0x20,0x20,0x1C,0x22,0x21 }, //173/ --> ­
  { 0x28,0x2C,0x2A,0x20,0x3F }, //174/ --> ®
  { 0x20,0x20,0x20,0x20,0x20 }, //175/ --> ¯
  { 0x20,0x20,0x20,0x1F,0x20 }, //176/ --> °
  { 0x30,0x40,0x40,0x3F,0x20 }, //177/ --> ±
  { 0x20,0x20,0x20,0x30,0x48 }, //178/ --> ²
  { 0x20,0x30,0x48,0x48,0x30 }, //179/ --> ³
  { 0x40,0x30,0x48,0x48,0x30 }, //180/ --> ´
  { 0x20,0x20,0x20,0x22,0x20 }, //181/ --> µ
  { 0x20,0x20,0x20,0x1A,0x20 }, //182/ --> ¶
  { 0x40,0x44,0x40,0x30,0x20 }, //183/ --> ·
  { 0x20,0x20,0x30,0x28,0x3A }, //184/ --> ¸
  { 0x20,0x30,0x28,0x3A,0x2C }, //185/ --> ¹
  { 0x18,0x14,0x14,0x18,0x20 }, //186/ --> º
  { 0x21,0x22,0x24,0x28,0x10 }, //187/ --> »
  { 0xB0,0xA8,0x78,0x20,0x20 }, //188/ --> ¼
  { 0x20,0x20,0xA0,0x20,0xA0 }, //189/ --> ½
  { 0x20,0x20,0xA0,0x30,0xA0 }, //190/ --> ¾
  { 0x60,0x80,0x80,0xA0,0x50 }, //191/ --> ¿
  { 0x1E,0x20,0x20,0x20,0x20 }, //192/ --> À
  { 0x20,0x30,0x28,0x28,0x20 }, //193/ --> Á
  { 0x04,0x02,0x02,0x3A,0x02 }, //194/ --> Â
  { 0x00,0x04,0x06,0x3D,0x05 }, //195/ --> Ã
  { 0x00,0x04,0xB6,0xAD,0x7D }, //196/ --> Ä
  { 0x00,0x80,0xC0,0xBF,0xA0 }, //197/ --> Å
  { 0x66,0x85,0x95,0xA8,0xA8 }, //198/ --> Æ
  { 0x00,0x00,0x3F,0x00,0x00 }, //199/ --> Ç
  { 0x30,0x28,0x20,0xA0,0x20 }, //200/ --> È
  { 0x00,0x30,0x2A,0x28,0x32 }, //201/ --> É
  { 0x30,0x28,0x22,0x20,0x22 }, //202/ --> Ê
  { 0x30,0x28,0x22,0x21,0x22 }, //203/ --> Ë
  { 0xC0,0xA8,0xA8,0x28,0xB0 }, //204/ --> Ì
  { 0xC0,0xA8,0xA8,0xA8,0x30 }, //205/ --> Í
  { 0xC0,0xA8,0xAA,0x28,0x30 }, //206/ --> Î
  { 0x00,0x24,0x24,0x24,0x38 }, //207/ --> Ï
  { 0x00,0x24,0x25,0x24,0x38 }, //208/ --> Ð
  { 0x80,0x80,0x40,0x30,0x00 }, //209/ --> Ñ
  { 0x00,0x80,0x80,0x40,0x34 }, //210/ --> Ò
  { 0x80,0x80,0x78,0x20,0x38 }, //211/ --> Ó
  { 0x80,0x80,0x78,0x22,0x39 }, //212/ --> Ô
  { 0x80,0x80,0x60,0x30,0x28 }, //213/ --> Õ
  { 0x80,0x80,0x60,0x30,0x28 }, //214/ --> Ö
  { 0x22,0x14,0x08,0x14,0x22 }, //215/ --> ×
  { 0x20,0x3E,0x30,0x28,0x28 }, //216/ --> Ø
  { 0x20,0x3E,0x30,0x28,0x2A }, //217/ --> Ù
  { 0x00,0x40,0xA0,0xB0,0x28 }, //218/ --> Ú
  { 0x00,0x40,0xA0,0xB0,0x2A }, //219/ --> Û
  { 0x20,0x20,0x20,0x20,0x20 }, //220/ --> Ü
  { 0x18,0x20,0x20,0x30,0x28 }, //221/ --> Ý
  { 0x60,0x80,0x80,0xB2,0xA8 }, //222/ --> Þ
  { 0x30,0x28,0x2C,0x2A,0x20 }, //223/ --> ß
  { 0x40,0xA9,0xAA,0xA8,0xF0 }, //224/ --> à
  { 0x00,0x60,0x80,0x80,0x7E }, //225/ --> á
  { 0x40,0xAA,0xA9,0xAA,0xF0 }, //226/ --> â
  { 0x00,0xC0,0x20,0x30,0x28 }, //227/ --> ã
  { 0x00,0x60,0x80,0x88,0x80 }, //228/ --> ä
  { 0x00,0x30,0x28,0x28,0x30 }, //229/ --> å
  { 0x00,0x00,0xB0,0xA8,0x78 }, //230/ --> æ
  { 0x26,0x25,0x25,0x28,0x10 }, //231/ --> ç
  { 0x22,0x22,0x26,0x29,0x10 }, //232/ --> è
  { 0x21,0x22,0x24,0xA8,0xD0 }, //233/ --> é
  { 0x70,0xAA,0xA9,0xAA,0x30 }, //234/ --> ê
  { 0x70,0xAA,0xA8,0xAA,0x30 }, //235/ --> ë
  { 0x30,0x40,0x40,0x50,0x28 }, //236/ --> ì
  { 0x30,0xC0,0x40,0xD0,0x28 }, //237/ --> í
  { 0x00,0x02,0x79,0x02,0x00 }, //238/ --> î
  { 0x00,0x02,0x78,0x02,0x00 }, //239/ --> ï
  { 0x00,0x00,0x00,0x00,0x05 }, //240/ --> ð
  { 0x00,0x00,0x04,0x03,0x0B }, //241/ --> ñ
  { 0xA0,0xA0,0x00,0x00,0x00 }, //242/ --> ò
  { 0x00,0x01,0x01,0x01,0x01 }, //243/ --> ó
  { 0x21,0x22,0x24,0x28,0x10 }, //244/ --> ô
  { 0x00,0x00,0x00,0x00,0x04 }, //245/ --> õ
  { 0x00,0x80,0x80,0x80,0x80 }, //246/ --> ö
  { 0x00,0x10,0x10,0x54,0x10 }, //247/ --> ÷
  { 0x00,0x02,0x04,0x02,0x04 }, //248/ --> ø
  { 0x24,0x26,0x25,0x25,0x20 }, //249/ --> ù
  { 0x24,0x26,0x25,0x35,0x20 }, //250/ --> ú
  { 0x6C,0x8A,0x8A,0xA0,0x50 }, //251/ --> û
  { 0xB6,0xAD,0x7D,0x24,0x20 }, //252/ --> ü
  { 0x19,0x14,0x15,0x18,0x20 }, //253/ --> ý
  { 0x02,0x02,0x1A,0x22,0x22 }, //254/ --> þ
0x00,0x40,0x60,0x50,0x48,0x50,0x40,0x40 //255/ --> ÿ
*/
};
