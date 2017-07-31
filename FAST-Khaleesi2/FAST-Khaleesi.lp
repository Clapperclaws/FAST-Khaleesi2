\ENCODING=ISO-8859-1
\Problem name: ilog.cplex

Minimize
 obj: w1144 + w1143 + w1142 + w1141 + w1140 + w1134 + w1133 + w1132 + w1131
      + w1130 + w1124 + w1123 + w1122 + w1121 + w1120 + w1114 + w1113 + w1112
      + w1111 + w1110 + w1104 + w1103 + w1102 + w1101 + w1100 + w1044 + w1043
      + w1042 + w1041 + w1040 + w1034 + w1033 + w1032 + w1031 + w1030 + w1024
      + w1023 + w1022 + w1021 + w1020 + w1014 + w1013 + w1012 + w1011 + w1010
      + w1004 + w1003 + w1002 + w1001 + w1000 + w944 + w943 + w942 + w941
      + w940 + w934 + w933 + w932 + w931 + w930 + w924 + w923 + w922 + w921
      + w920 + w914 + w913 + w912 + w911 + w910 + w904 + w903 + w902 + w901
      + w900 + w844 + w843 + w842 + w841 + w840 + w834 + w833 + w832 + w831
      + w830 + w824 + w823 + w822 + w821 + w820 + w814 + w813 + w812 + w811
      + w810 + w804 + w803 + w802 + w801 + w800 + w744 + w743 + w742 + w741
      + w740 + w734 + w733 + w732 + w731 + w730 + w724 + w723 + w722 + w721
      + w720 + w714 + w713 + w712 + w711 + w710 + w704 + w703 + w702 + w701
      + w700 + w644 + w643 + w642 + w641 + w640 + w634 + w633 + w632 + w631
      + w630 + w624 + w623 + w622 + w621 + w620 + w614 + w613 + w612 + w611
      + w610 + w604 + w603 + w602 + w601 + w600 + w544 + w543 + w542 + w541
      + w540 + w534 + w533 + w532 + w531 + w530 + w524 + w523 + w522 + w521
      + w520 + w514 + w513 + w512 + w511 + w510 + w504 + w503 + w502 + w501
      + w500 + w444 + w443 + w442 + w441 + w440 + w434 + w433 + w432 + w431
      + w430 + w424 + w423 + w422 + w421 + w420 + w414 + w413 + w412 + w411
      + w410 + w404 + w403 + w402 + w401 + w400 + w344 + w343 + w342 + w341
      + w340 + w334 + w333 + w332 + w331 + w330 + w324 + w323 + w322 + w321
      + w320 + w314 + w313 + w312 + w311 + w310 + w304 + w303 + w302 + w301
      + w300 + w244 + w243 + w242 + w241 + w240 + w234 + w233 + w232 + w231
      + w230 + w224 + w223 + w222 + w221 + w220 + w214 + w213 + w212 + w211
      + w210 + w204 + w203 + w202 + w201 + w200 + w144 + w143 + w142 + w141
      + w140 + w134 + w133 + w132 + w131 + w130 + w124 + w123 + w122 + w121
      + w120 + w114 + w113 + w112 + w111 + w110 + w104 + w103 + w102 + w101
      + w100 + w044 + w043 + w042 + w041 + w040 + w034 + w033 + w032 + w031
      + w030 + w024 + w023 + w022 + w021 + w020 + w014 + w013 + w012 + w011
      + w010 + w004 + w003 + w002 + w001 + w000
Subject To
 NF_0PlacementConstraint#0:      theta04 + theta03 + theta02 + theta01
                                 + theta00  = 1
 NF_1PlacementConstraint#1:      theta14 + theta13 + theta12 + theta11
                                 + theta10  = 1
 NF_2PlacementConstraint#2:      theta24 + theta23 + theta22 + theta21
                                 + theta20  = 1
 NF_3PlacementConstraint#3:      theta34 + theta33 + theta32 + theta31
                                 + theta30  = 1
 Source_of_E_0_Placement_Constraint#4: - theta00 + x00 <= 0
 Destination_of_E_0_Placement_Constraint#5: - theta10 + y00 <= 0
 Source_of_E_0_Placement_Constraint#6: - theta01 + x01 <= 0
 Destination_of_E_0_Placement_Constraint#7: - theta11 + y01 <= 0
 Source_of_E_0_Placement_Constraint#8: - theta02 + x02 <= 0
 Destination_of_E_0_Placement_Constraint#9: - theta12 + y02 <= 0
 Source_of_E_0_Placement_Constraint#10: - theta03 + x03 <= 0
 Destination_of_E_0_Placement_Constraint#11: - theta13 + y03 <= 0
 Source_of_E_0_Placement_Constraint#12: - theta04 + x04 <= 0
 Destination_of_E_0_Placement_Constraint#13: - theta14 + y04 <= 0
 Source_of_E_1_Placement_Constraint#14: - theta10 + x10 <= 0
 Destination_of_E_1_Placement_Constraint#15: - theta20 + y10 <= 0
 Source_of_E_1_Placement_Constraint#16: - theta11 + x11 <= 0
 Destination_of_E_1_Placement_Constraint#17: - theta21 + y11 <= 0
 Source_of_E_1_Placement_Constraint#18: - theta12 + x12 <= 0
 Destination_of_E_1_Placement_Constraint#19: - theta22 + y12 <= 0
 Source_of_E_1_Placement_Constraint#20: - theta13 + x13 <= 0
 Destination_of_E_1_Placement_Constraint#21: - theta23 + y13 <= 0
 Source_of_E_1_Placement_Constraint#22: - theta14 + x14 <= 0
 Destination_of_E_1_Placement_Constraint#23: - theta24 + y14 <= 0
 Source_of_E_2_Placement_Constraint#24: - theta20 + x20 <= 0
 Destination_of_E_2_Placement_Constraint#25: - theta30 + y20 <= 0
 Source_of_E_2_Placement_Constraint#26: - theta21 + x21 <= 0
 Destination_of_E_2_Placement_Constraint#27: - theta31 + y21 <= 0
 Source_of_E_2_Placement_Constraint#28: - theta22 + x22 <= 0
 Destination_of_E_2_Placement_Constraint#29: - theta32 + y22 <= 0
 Source_of_E_2_Placement_Constraint#30: - theta23 + x23 <= 0
 Destination_of_E_2_Placement_Constraint#31: - theta33 + y23 <= 0
 Source_of_E_2_Placement_Constraint#32: - theta24 + x24 <= 0
 Destination_of_E_2_Placement_Constraint#33: - theta34 + y24 <= 0
 Source_of_E_3_Placement_Constraint#34: - theta10 + x30 <= 0
 Destination_of_E_3_Placement_Constraint#35: - theta00 + y30 <= 0
 Source_of_E_3_Placement_Constraint#36: - theta11 + x31 <= 0
 Destination_of_E_3_Placement_Constraint#37: - theta01 + y31 <= 0
 Source_of_E_3_Placement_Constraint#38: - theta12 + x32 <= 0
 Destination_of_E_3_Placement_Constraint#39: - theta02 + y32 <= 0
 Source_of_E_3_Placement_Constraint#40: - theta13 + x33 <= 0
 Destination_of_E_3_Placement_Constraint#41: - theta03 + y33 <= 0
 Source_of_E_3_Placement_Constraint#42: - theta14 + x34 <= 0
 Destination_of_E_3_Placement_Constraint#43: - theta04 + y34 <= 0
 Source_of_E_4_Placement_Constraint#44: - theta00 + x40 <= 0
 Destination_of_E_4_Placement_Constraint#45: - theta20 + y40 <= 0
 Source_of_E_4_Placement_Constraint#46: - theta01 + x41 <= 0
 Destination_of_E_4_Placement_Constraint#47: - theta21 + y41 <= 0
 Source_of_E_4_Placement_Constraint#48: - theta02 + x42 <= 0
 Destination_of_E_4_Placement_Constraint#49: - theta22 + y42 <= 0
 Source_of_E_4_Placement_Constraint#50: - theta03 + x43 <= 0
 Destination_of_E_4_Placement_Constraint#51: - theta23 + y43 <= 0
 Source_of_E_4_Placement_Constraint#52: - theta04 + x44 <= 0
 Destination_of_E_4_Placement_Constraint#53: - theta24 + y44 <= 0
 Source_of_E_5_Placement_Constraint#54: - theta20 + x50 <= 0
 Destination_of_E_5_Placement_Constraint#55: - theta00 + y50 <= 0
 Source_of_E_5_Placement_Constraint#56: - theta21 + x51 <= 0
 Destination_of_E_5_Placement_Constraint#57: - theta01 + y51 <= 0
 Source_of_E_5_Placement_Constraint#58: - theta22 + x52 <= 0
 Destination_of_E_5_Placement_Constraint#59: - theta02 + y52 <= 0
 Source_of_E_5_Placement_Constraint#60: - theta23 + x53 <= 0
 Destination_of_E_5_Placement_Constraint#61: - theta03 + y53 <= 0
 Source_of_E_5_Placement_Constraint#62: - theta24 + x54 <= 0
 Destination_of_E_5_Placement_Constraint#63: - theta04 + y54 <= 0
 Source_of_E_6_Placement_Constraint#64: - theta00 + x60 <= 0
 Destination_of_E_6_Placement_Constraint#65: - theta30 + y60 <= 0
 Source_of_E_6_Placement_Constraint#66: - theta01 + x61 <= 0
 Destination_of_E_6_Placement_Constraint#67: - theta31 + y61 <= 0
 Source_of_E_6_Placement_Constraint#68: - theta02 + x62 <= 0
 Destination_of_E_6_Placement_Constraint#69: - theta32 + y62 <= 0
 Source_of_E_6_Placement_Constraint#70: - theta03 + x63 <= 0
 Destination_of_E_6_Placement_Constraint#71: - theta33 + y63 <= 0
 Source_of_E_6_Placement_Constraint#72: - theta04 + x64 <= 0
 Destination_of_E_6_Placement_Constraint#73: - theta34 + y64 <= 0
 Source_of_E_7_Placement_Constraint#74: - theta30 + x70 <= 0
 Destination_of_E_7_Placement_Constraint#75: - theta00 + y70 <= 0
 Source_of_E_7_Placement_Constraint#76: - theta31 + x71 <= 0
 Destination_of_E_7_Placement_Constraint#77: - theta01 + y71 <= 0
 Source_of_E_7_Placement_Constraint#78: - theta32 + x72 <= 0
 Destination_of_E_7_Placement_Constraint#79: - theta02 + y72 <= 0
 Source_of_E_7_Placement_Constraint#80: - theta33 + x73 <= 0
 Destination_of_E_7_Placement_Constraint#81: - theta03 + y73 <= 0
 Source_of_E_7_Placement_Constraint#82: - theta34 + x74 <= 0
 Destination_of_E_7_Placement_Constraint#83: - theta04 + y74 <= 0
 Source_of_E_8_Placement_Constraint#84: - theta20 + x80 <= 0
 Destination_of_E_8_Placement_Constraint#85: - theta10 + y80 <= 0
 Source_of_E_8_Placement_Constraint#86: - theta21 + x81 <= 0
 Destination_of_E_8_Placement_Constraint#87: - theta11 + y81 <= 0
 Source_of_E_8_Placement_Constraint#88: - theta22 + x82 <= 0
 Destination_of_E_8_Placement_Constraint#89: - theta12 + y82 <= 0
 Source_of_E_8_Placement_Constraint#90: - theta23 + x83 <= 0
 Destination_of_E_8_Placement_Constraint#91: - theta13 + y83 <= 0
 Source_of_E_8_Placement_Constraint#92: - theta24 + x84 <= 0
 Destination_of_E_8_Placement_Constraint#93: - theta14 + y84 <= 0
 Source_of_E_9_Placement_Constraint#94: - theta10 + x90 <= 0
 Destination_of_E_9_Placement_Constraint#95: - theta30 + y90 <= 0
 Source_of_E_9_Placement_Constraint#96: - theta11 + x91 <= 0
 Destination_of_E_9_Placement_Constraint#97: - theta31 + y91 <= 0
 Source_of_E_9_Placement_Constraint#98: - theta12 + x92 <= 0
 Destination_of_E_9_Placement_Constraint#99: - theta32 + y92 <= 0
 Source_of_E_9_Placement_Constraint#100: - theta13 + x93 <= 0
 Destination_of_E_9_Placement_Constraint#101: - theta33 + y93 <= 0
 Source_of_E_9_Placement_Constraint#102: - theta14 + x94 <= 0
 Destination_of_E_9_Placement_Constraint#103: - theta34 + y94 <= 0
 Source_of_E_10_Placement_Constraint#104: - theta30 + x100 <= 0
 Destination_of_E_10_Placement_Constraint#105: - theta10 + y100 <= 0
 Source_of_E_10_Placement_Constraint#106: - theta31 + x101 <= 0
 Destination_of_E_10_Placement_Constraint#107: - theta11 + y101 <= 0
 Source_of_E_10_Placement_Constraint#108: - theta32 + x102 <= 0
 Destination_of_E_10_Placement_Constraint#109: - theta12 + y102 <= 0
 Source_of_E_10_Placement_Constraint#110: - theta33 + x103 <= 0
 Destination_of_E_10_Placement_Constraint#111: - theta13 + y103 <= 0
 Source_of_E_10_Placement_Constraint#112: - theta34 + x104 <= 0
 Destination_of_E_10_Placement_Constraint#113: - theta14 + y104 <= 0
 Source_of_E_11_Placement_Constraint#114: - theta30 + x110 <= 0
 Destination_of_E_11_Placement_Constraint#115: - theta20 + y110 <= 0
 Source_of_E_11_Placement_Constraint#116: - theta31 + x111 <= 0
 Destination_of_E_11_Placement_Constraint#117: - theta21 + y111 <= 0
 Source_of_E_11_Placement_Constraint#118: - theta32 + x112 <= 0
 Destination_of_E_11_Placement_Constraint#119: - theta22 + y112 <= 0
 Source_of_E_11_Placement_Constraint#120: - theta33 + x113 <= 0
 Destination_of_E_11_Placement_Constraint#121: - theta23 + y113 <= 0
 Source_of_E_11_Placement_Constraint#122: - theta34 + x114 <= 0
 Destination_of_E_11_Placement_Constraint#123: - theta24 + y114 <= 0
 CapacityConstraints_Node0#124:  theta00 + theta10 + theta20 + theta30 <= 1
 CapacityConstraints_Node1#125:  theta01 + theta11 + theta21 + theta31 <= 1
 CapacityConstraints_Node2#126:  theta02 + theta12 + theta22 + theta32 <= 1
 CapacityConstraints_Node3#127:  theta03 + theta13 + theta23 + theta33 <= 1
 CapacityConstraints_Node4#128:  theta04 + theta14 + theta24 + theta34 <= 1
 LinkPlacementConstraint0#129:   - 0.5 x00 - 0.5 y00 - 0.5 x01 - 0.5 y01
                                 - 0.5 x02 - 0.5 y02 - 0.5 x03 - 0.5 y03
                                 - 0.5 x04 - 0.5 y04 + z0  = 0
 LinkPlacementConstraint1#130:   - 0.5 x10 - 0.5 y10 - 0.5 x11 - 0.5 y11
                                 - 0.5 x12 - 0.5 y12 - 0.5 x13 - 0.5 y13
                                 - 0.5 x14 - 0.5 y14 + z1  = 0
 LinkPlacementConstraint2#131:   - 0.5 x20 - 0.5 y20 - 0.5 x21 - 0.5 y21
                                 - 0.5 x22 - 0.5 y22 - 0.5 x23 - 0.5 y23
                                 - 0.5 x24 - 0.5 y24 + z2  = 0
 LinkPlacementConstraint3#132:   - 0.5 x30 - 0.5 y30 - 0.5 x31 - 0.5 y31
                                 - 0.5 x32 - 0.5 y32 - 0.5 x33 - 0.5 y33
                                 - 0.5 x34 - 0.5 y34 + z3  = 0
 LinkPlacementConstraint4#133:   - 0.5 x40 - 0.5 y40 - 0.5 x41 - 0.5 y41
                                 - 0.5 x42 - 0.5 y42 - 0.5 x43 - 0.5 y43
                                 - 0.5 x44 - 0.5 y44 + z4  = 0
 LinkPlacementConstraint5#134:   - 0.5 x50 - 0.5 y50 - 0.5 x51 - 0.5 y51
                                 - 0.5 x52 - 0.5 y52 - 0.5 x53 - 0.5 y53
                                 - 0.5 x54 - 0.5 y54 + z5  = 0
 LinkPlacementConstraint6#135:   - 0.5 x60 - 0.5 y60 - 0.5 x61 - 0.5 y61
                                 - 0.5 x62 - 0.5 y62 - 0.5 x63 - 0.5 y63
                                 - 0.5 x64 - 0.5 y64 + z6  = 0
 LinkPlacementConstraint7#136:   - 0.5 x70 - 0.5 y70 - 0.5 x71 - 0.5 y71
                                 - 0.5 x72 - 0.5 y72 - 0.5 x73 - 0.5 y73
                                 - 0.5 x74 - 0.5 y74 + z7  = 0
 LinkPlacementConstraint8#137:   - 0.5 x80 - 0.5 y80 - 0.5 x81 - 0.5 y81
                                 - 0.5 x82 - 0.5 y82 - 0.5 x83 - 0.5 y83
                                 - 0.5 x84 - 0.5 y84 + z8  = 0
 LinkPlacementConstraint9#138:   - 0.5 x90 - 0.5 y90 - 0.5 x91 - 0.5 y91
                                 - 0.5 x92 - 0.5 y92 - 0.5 x93 - 0.5 y93
                                 - 0.5 x94 - 0.5 y94 + z9  = 0
 LinkPlacementConstraint10#139:  - 0.5 x100 - 0.5 y100 - 0.5 x101 - 0.5 y101
                                 - 0.5 x102 - 0.5 y102 - 0.5 x103 - 0.5 y103
                                 - 0.5 x104 - 0.5 y104 + z10  = 0
 LinkPlacementConstraint11#140:  - 0.5 x110 - 0.5 y110 - 0.5 x111 - 0.5 y111
                                 - 0.5 x112 - 0.5 y112 - 0.5 x113 - 0.5 y113
                                 - 0.5 x114 - 0.5 y114 + z11  = 0
 IloC141#141:                    z0 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8
                                 + z9 + z10 + z11  = 3
 IloC142#142:                    w010 + w001 <= 1
 IloC143#143:                    w020 + w002 <= 1
 Routing_Constraint_0#144:       - w020 - w010 + w002 + w001 - x00 + y00  = 0
 IloC145#145:                    w010 + w001 <= 1
 IloC146#146:                    w021 + w012 <= 1
 Routing_Constraint_0#147:       - w021 + w012 + w010 - w001 - x01 + y01  = 0
 IloC148#148:                    w020 + w002 <= 1
 IloC149#149:                    w021 + w012 <= 1
 IloC150#150:                    w032 + w023 <= 1
 Routing_Constraint_0#151:       - w032 + w023 + w021 + w020 - w012 - w002
                                 - x02 + y02  = 0
 IloC152#152:                    w032 + w023 <= 1
 IloC153#153:                    w043 + w034 <= 1
 Routing_Constraint_0#154:       - w043 + w034 + w032 - w023 - x03 + y03  = 0
 IloC155#155:                    w043 + w034 <= 1
 Routing_Constraint_0#156:       w043 - w034 - x04 + y04  = 0
 IloC157#157:                    w110 + w101 <= 1
 IloC158#158:                    w120 + w102 <= 1
 Routing_Constraint_1#159:       - w120 - w110 + w102 + w101 - x10 + y10  = 0
 IloC160#160:                    w110 + w101 <= 1
 IloC161#161:                    w121 + w112 <= 1
 Routing_Constraint_1#162:       - w121 + w112 + w110 - w101 - x11 + y11  = 0
 IloC163#163:                    w120 + w102 <= 1
 IloC164#164:                    w121 + w112 <= 1
 IloC165#165:                    w132 + w123 <= 1
 Routing_Constraint_1#166:       - w132 + w123 + w121 + w120 - w112 - w102
                                 - x12 + y12  = 0
 IloC167#167:                    w132 + w123 <= 1
 IloC168#168:                    w143 + w134 <= 1
 Routing_Constraint_1#169:       - w143 + w134 + w132 - w123 - x13 + y13  = 0
 IloC170#170:                    w143 + w134 <= 1
 Routing_Constraint_1#171:       w143 - w134 - x14 + y14  = 0
 IloC172#172:                    w210 + w201 <= 1
 IloC173#173:                    w220 + w202 <= 1
 Routing_Constraint_2#174:       - w220 - w210 + w202 + w201 - x20 + y20  = 0
 IloC175#175:                    w210 + w201 <= 1
 IloC176#176:                    w221 + w212 <= 1
 Routing_Constraint_2#177:       - w221 + w212 + w210 - w201 - x21 + y21  = 0
 IloC178#178:                    w220 + w202 <= 1
 IloC179#179:                    w221 + w212 <= 1
 IloC180#180:                    w232 + w223 <= 1
 Routing_Constraint_2#181:       - w232 + w223 + w221 + w220 - w212 - w202
                                 - x22 + y22  = 0
 IloC182#182:                    w232 + w223 <= 1
 IloC183#183:                    w243 + w234 <= 1
 Routing_Constraint_2#184:       - w243 + w234 + w232 - w223 - x23 + y23  = 0
 IloC185#185:                    w243 + w234 <= 1
 Routing_Constraint_2#186:       w243 - w234 - x24 + y24  = 0
 IloC187#187:                    w310 + w301 <= 1
 IloC188#188:                    w320 + w302 <= 1
 Routing_Constraint_3#189:       - w320 - w310 + w302 + w301 - x30 + y30  = 0
 IloC190#190:                    w310 + w301 <= 1
 IloC191#191:                    w321 + w312 <= 1
 Routing_Constraint_3#192:       - w321 + w312 + w310 - w301 - x31 + y31  = 0
 IloC193#193:                    w320 + w302 <= 1
 IloC194#194:                    w321 + w312 <= 1
 IloC195#195:                    w332 + w323 <= 1
 Routing_Constraint_3#196:       - w332 + w323 + w321 + w320 - w312 - w302
                                 - x32 + y32  = 0
 IloC197#197:                    w332 + w323 <= 1
 IloC198#198:                    w343 + w334 <= 1
 Routing_Constraint_3#199:       - w343 + w334 + w332 - w323 - x33 + y33  = 0
 IloC200#200:                    w343 + w334 <= 1
 Routing_Constraint_3#201:       w343 - w334 - x34 + y34  = 0
 IloC202#202:                    w410 + w401 <= 1
 IloC203#203:                    w420 + w402 <= 1
 Routing_Constraint_4#204:       - w420 - w410 + w402 + w401 - x40 + y40  = 0
 IloC205#205:                    w410 + w401 <= 1
 IloC206#206:                    w421 + w412 <= 1
 Routing_Constraint_4#207:       - w421 + w412 + w410 - w401 - x41 + y41  = 0
 IloC208#208:                    w420 + w402 <= 1
 IloC209#209:                    w421 + w412 <= 1
 IloC210#210:                    w432 + w423 <= 1
 Routing_Constraint_4#211:       - w432 + w423 + w421 + w420 - w412 - w402
                                 - x42 + y42  = 0
 IloC212#212:                    w432 + w423 <= 1
 IloC213#213:                    w443 + w434 <= 1
 Routing_Constraint_4#214:       - w443 + w434 + w432 - w423 - x43 + y43  = 0
 IloC215#215:                    w443 + w434 <= 1
 Routing_Constraint_4#216:       w443 - w434 - x44 + y44  = 0
 IloC217#217:                    w510 + w501 <= 1
 IloC218#218:                    w520 + w502 <= 1
 Routing_Constraint_5#219:       - w520 - w510 + w502 + w501 - x50 + y50  = 0
 IloC220#220:                    w510 + w501 <= 1
 IloC221#221:                    w521 + w512 <= 1
 Routing_Constraint_5#222:       - w521 + w512 + w510 - w501 - x51 + y51  = 0
 IloC223#223:                    w520 + w502 <= 1
 IloC224#224:                    w521 + w512 <= 1
 IloC225#225:                    w532 + w523 <= 1
 Routing_Constraint_5#226:       - w532 + w523 + w521 + w520 - w512 - w502
                                 - x52 + y52  = 0
 IloC227#227:                    w532 + w523 <= 1
 IloC228#228:                    w543 + w534 <= 1
 Routing_Constraint_5#229:       - w543 + w534 + w532 - w523 - x53 + y53  = 0
 IloC230#230:                    w543 + w534 <= 1
 Routing_Constraint_5#231:       w543 - w534 - x54 + y54  = 0
 IloC232#232:                    w610 + w601 <= 1
 IloC233#233:                    w620 + w602 <= 1
 Routing_Constraint_6#234:       - w620 - w610 + w602 + w601 - x60 + y60  = 0
 IloC235#235:                    w610 + w601 <= 1
 IloC236#236:                    w621 + w612 <= 1
 Routing_Constraint_6#237:       - w621 + w612 + w610 - w601 - x61 + y61  = 0
 IloC238#238:                    w620 + w602 <= 1
 IloC239#239:                    w621 + w612 <= 1
 IloC240#240:                    w632 + w623 <= 1
 Routing_Constraint_6#241:       - w632 + w623 + w621 + w620 - w612 - w602
                                 - x62 + y62  = 0
 IloC242#242:                    w632 + w623 <= 1
 IloC243#243:                    w643 + w634 <= 1
 Routing_Constraint_6#244:       - w643 + w634 + w632 - w623 - x63 + y63  = 0
 IloC245#245:                    w643 + w634 <= 1
 Routing_Constraint_6#246:       w643 - w634 - x64 + y64  = 0
 IloC247#247:                    w710 + w701 <= 1
 IloC248#248:                    w720 + w702 <= 1
 Routing_Constraint_7#249:       - w720 - w710 + w702 + w701 - x70 + y70  = 0
 IloC250#250:                    w710 + w701 <= 1
 IloC251#251:                    w721 + w712 <= 1
 Routing_Constraint_7#252:       - w721 + w712 + w710 - w701 - x71 + y71  = 0
 IloC253#253:                    w720 + w702 <= 1
 IloC254#254:                    w721 + w712 <= 1
 IloC255#255:                    w732 + w723 <= 1
 Routing_Constraint_7#256:       - w732 + w723 + w721 + w720 - w712 - w702
                                 - x72 + y72  = 0
 IloC257#257:                    w732 + w723 <= 1
 IloC258#258:                    w743 + w734 <= 1
 Routing_Constraint_7#259:       - w743 + w734 + w732 - w723 - x73 + y73  = 0
 IloC260#260:                    w743 + w734 <= 1
 Routing_Constraint_7#261:       w743 - w734 - x74 + y74  = 0
 IloC262#262:                    w810 + w801 <= 1
 IloC263#263:                    w820 + w802 <= 1
 Routing_Constraint_8#264:       - w820 - w810 + w802 + w801 - x80 + y80  = 0
 IloC265#265:                    w810 + w801 <= 1
 IloC266#266:                    w821 + w812 <= 1
 Routing_Constraint_8#267:       - w821 + w812 + w810 - w801 - x81 + y81  = 0
 IloC268#268:                    w820 + w802 <= 1
 IloC269#269:                    w821 + w812 <= 1
 IloC270#270:                    w832 + w823 <= 1
 Routing_Constraint_8#271:       - w832 + w823 + w821 + w820 - w812 - w802
                                 - x82 + y82  = 0
 IloC272#272:                    w832 + w823 <= 1
 IloC273#273:                    w843 + w834 <= 1
 Routing_Constraint_8#274:       - w843 + w834 + w832 - w823 - x83 + y83  = 0
 IloC275#275:                    w843 + w834 <= 1
 Routing_Constraint_8#276:       w843 - w834 - x84 + y84  = 0
 IloC277#277:                    w910 + w901 <= 1
 IloC278#278:                    w920 + w902 <= 1
 Routing_Constraint_9#279:       - w920 - w910 + w902 + w901 - x90 + y90  = 0
 IloC280#280:                    w910 + w901 <= 1
 IloC281#281:                    w921 + w912 <= 1
 Routing_Constraint_9#282:       - w921 + w912 + w910 - w901 - x91 + y91  = 0
 IloC283#283:                    w920 + w902 <= 1
 IloC284#284:                    w921 + w912 <= 1
 IloC285#285:                    w932 + w923 <= 1
 Routing_Constraint_9#286:       - w932 + w923 + w921 + w920 - w912 - w902
                                 - x92 + y92  = 0
 IloC287#287:                    w932 + w923 <= 1
 IloC288#288:                    w943 + w934 <= 1
 Routing_Constraint_9#289:       - w943 + w934 + w932 - w923 - x93 + y93  = 0
 IloC290#290:                    w943 + w934 <= 1
 Routing_Constraint_9#291:       w943 - w934 - x94 + y94  = 0
 IloC292#292:                    w1010 + w1001 <= 1
 IloC293#293:                    w1020 + w1002 <= 1
 Routing_Constraint_10#294:      - w1020 - w1010 + w1002 + w1001 - x100 + y100
                                  = 0
 IloC295#295:                    w1010 + w1001 <= 1
 IloC296#296:                    w1021 + w1012 <= 1
 Routing_Constraint_10#297:      - w1021 + w1012 + w1010 - w1001 - x101 + y101
                                  = 0
 IloC298#298:                    w1020 + w1002 <= 1
 IloC299#299:                    w1021 + w1012 <= 1
 IloC300#300:                    w1032 + w1023 <= 1
 Routing_Constraint_10#301:      - w1032 + w1023 + w1021 + w1020 - w1012
                                 - w1002 - x102 + y102  = 0
 IloC302#302:                    w1032 + w1023 <= 1
 IloC303#303:                    w1043 + w1034 <= 1
 Routing_Constraint_10#304:      - w1043 + w1034 + w1032 - w1023 - x103 + y103
                                  = 0
 IloC305#305:                    w1043 + w1034 <= 1
 Routing_Constraint_10#306:      w1043 - w1034 - x104 + y104  = 0
 IloC307#307:                    w1110 + w1101 <= 1
 IloC308#308:                    w1120 + w1102 <= 1
 Routing_Constraint_11#309:      - w1120 - w1110 + w1102 + w1101 - x110 + y110
                                  = 0
 IloC310#310:                    w1110 + w1101 <= 1
 IloC311#311:                    w1121 + w1112 <= 1
 Routing_Constraint_11#312:      - w1121 + w1112 + w1110 - w1101 - x111 + y111
                                  = 0
 IloC313#313:                    w1120 + w1102 <= 1
 IloC314#314:                    w1121 + w1112 <= 1
 IloC315#315:                    w1132 + w1123 <= 1
 Routing_Constraint_11#316:      - w1132 + w1123 + w1121 + w1120 - w1112
                                 - w1102 - x112 + y112  = 0
 IloC317#317:                    w1132 + w1123 <= 1
 IloC318#318:                    w1143 + w1134 <= 1
 Routing_Constraint_11#319:      - w1143 + w1134 + w1132 - w1123 - x113 + y113
                                  = 0
 IloC320#320:                    w1143 + w1134 <= 1
 Routing_Constraint_11#321:      w1143 - w1134 - x114 + y114  = 0
 A_single_VLink_originating_from_NF0#322: x00 + x01 + x02 + x03 + x04 + x40
                                 + x41 + x42 + x43 + x44 + x60 + x61 + x62
                                 + x63 + x64 <= 1
 A_single_VLink_destined_to_NF0#323: y30 + y31 + y32 + y33 + y34 + y50 + y51
                                 + y52 + y53 + y54 + y70 + y71 + y72 + y73
                                 + y74 <= 1
 A_single_VLink_originating_from_NF1#324: x10 + x11 + x12 + x13 + x14 + x30
                                 + x31 + x32 + x33 + x34 + x90 + x91 + x92
                                 + x93 + x94 <= 1
 A_single_VLink_destined_to_NF1#325: y00 + y01 + y02 + y03 + y04 + y80 + y81
                                 + y82 + y83 + y84 + y100 + y101 + y102 + y103
                                 + y104 <= 1
 A_single_VLink_originating_from_NF2#326: x20 + x21 + x22 + x23 + x24 + x50
                                 + x51 + x52 + x53 + x54 + x80 + x81 + x82
                                 + x83 + x84 <= 1
 A_single_VLink_destined_to_NF2#327: y10 + y11 + y12 + y13 + y14 + y40 + y41
                                 + y42 + y43 + y44 + y110 + y111 + y112 + y113
                                 + y114 <= 1
 A_single_VLink_originating_from_NF3#328: x70 + x71 + x72 + x73 + x74 + x100
                                 + x101 + x102 + x103 + x104 + x110 + x111
                                 + x112 + x113 + x114 <= 1
 A_single_VLink_destined_to_NF3#329: y20 + y21 + y22 + y23 + y24 + y60 + y61
                                 + y62 + y63 + y64 + y90 + y91 + y92 + y93
                                 + y94 <= 1
 A_single_VLink_originating_from_NF0#330: z0 + z4 + z6 <= 1
 A_single_VLink_destined_to_NF0#331: z3 + z5 + z7 <= 1
 A_single_VLink_originating_from_NF1#332: z1 + z3 + z9 <= 1
 A_single_VLink_destined_to_NF1#333: z0 + z8 + z10 <= 1
 A_single_VLink_originating_from_NF2#334: z2 + z5 + z8 <= 1
 A_single_VLink_destined_to_NF2#335: z1 + z4 + z11 <= 1
 A_single_VLink_originating_from_NF3#336: z7 + z10 + z11 <= 1
 A_single_VLink_destined_to_NF3#337: z2 + z6 + z9 <= 1
 LinkCapacity01#338:             w1110 + w1101 + w1010 + w1001 + w910 + w901
                                 + w810 + w801 + w710 + w701 + w610 + w601
                                 + w510 + w501 + w410 + w401 + w310 + w301
                                 + w210 + w201 + w110 + w101 + w010 + w001
                                 <= 10
 LinkCapacity02#339:             w1120 + w1102 + w1020 + w1002 + w920 + w902
                                 + w820 + w802 + w720 + w702 + w620 + w602
                                 + w520 + w502 + w420 + w402 + w320 + w302
                                 + w220 + w202 + w120 + w102 + w020 + w002
                                 <= 10
 LinkCapacity12#340:             w1121 + w1112 + w1021 + w1012 + w921 + w912
                                 + w821 + w812 + w721 + w712 + w621 + w612
                                 + w521 + w512 + w421 + w412 + w321 + w312
                                 + w221 + w212 + w121 + w112 + w021 + w012
                                 <= 10
 LinkCapacity23#341:             w1132 + w1123 + w1032 + w1023 + w932 + w923
                                 + w832 + w823 + w732 + w723 + w632 + w623
                                 + w532 + w523 + w432 + w423 + w332 + w323
                                 + w232 + w223 + w132 + w123 + w032 + w023
                                 <= 10
 LinkCapacity34#342:             w1143 + w1134 + w1043 + w1034 + w943 + w934
                                 + w843 + w834 + w743 + w734 + w643 + w634
                                 + w543 + w534 + w443 + w434 + w343 + w334
                                 + w243 + w234 + w143 + w134 + w043 + w034
                                 <= 10
 chain_semantic03#343:           z0 + z3 <= 1
 chain_semantic18#344:           z1 + z8 <= 1
 chain_semantic211#345:          z2 + z11 <= 1
 chain_semantic30#346:           z0 + z3 <= 1
 chain_semantic45#347:           z4 + z5 <= 1
 chain_semantic54#348:           z4 + z5 <= 1
 chain_semantic67#349:           z6 + z7 <= 1
 chain_semantic76#350:           z6 + z7 <= 1
 chain_semantic81#351:           z1 + z8 <= 1
 chain_semantic910#352:          z9 + z10 <= 1
 chain_semantic109#353:          z9 + z10 <= 1
 chain_semantic112#354:          z2 + z11 <= 1
Bounds
 0 <= w1144 <= 1
 0 <= w1143 <= 1
 0 <= w1142 <= 1
 0 <= w1141 <= 1
 0 <= w1140 <= 1
 0 <= w1134 <= 1
 0 <= w1133 <= 1
 0 <= w1132 <= 1
 0 <= w1131 <= 1
 0 <= w1130 <= 1
 0 <= w1124 <= 1
 0 <= w1123 <= 1
 0 <= w1122 <= 1
 0 <= w1121 <= 1
 0 <= w1120 <= 1
 0 <= w1114 <= 1
 0 <= w1113 <= 1
 0 <= w1112 <= 1
 0 <= w1111 <= 1
 0 <= w1110 <= 1
 0 <= w1104 <= 1
 0 <= w1103 <= 1
 0 <= w1102 <= 1
 0 <= w1101 <= 1
 0 <= w1100 <= 1
 0 <= w1044 <= 1
 0 <= w1043 <= 1
 0 <= w1042 <= 1
 0 <= w1041 <= 1
 0 <= w1040 <= 1
 0 <= w1034 <= 1
 0 <= w1033 <= 1
 0 <= w1032 <= 1
 0 <= w1031 <= 1
 0 <= w1030 <= 1
 0 <= w1024 <= 1
 0 <= w1023 <= 1
 0 <= w1022 <= 1
 0 <= w1021 <= 1
 0 <= w1020 <= 1
 0 <= w1014 <= 1
 0 <= w1013 <= 1
 0 <= w1012 <= 1
 0 <= w1011 <= 1
 0 <= w1010 <= 1
 0 <= w1004 <= 1
 0 <= w1003 <= 1
 0 <= w1002 <= 1
 0 <= w1001 <= 1
 0 <= w1000 <= 1
 0 <= w944 <= 1
 0 <= w943 <= 1
 0 <= w942 <= 1
 0 <= w941 <= 1
 0 <= w940 <= 1
 0 <= w934 <= 1
 0 <= w933 <= 1
 0 <= w932 <= 1
 0 <= w931 <= 1
 0 <= w930 <= 1
 0 <= w924 <= 1
 0 <= w923 <= 1
 0 <= w922 <= 1
 0 <= w921 <= 1
 0 <= w920 <= 1
 0 <= w914 <= 1
 0 <= w913 <= 1
 0 <= w912 <= 1
 0 <= w911 <= 1
 0 <= w910 <= 1
 0 <= w904 <= 1
 0 <= w903 <= 1
 0 <= w902 <= 1
 0 <= w901 <= 1
 0 <= w900 <= 1
 0 <= w844 <= 1
 0 <= w843 <= 1
 0 <= w842 <= 1
 0 <= w841 <= 1
 0 <= w840 <= 1
 0 <= w834 <= 1
 0 <= w833 <= 1
 0 <= w832 <= 1
 0 <= w831 <= 1
 0 <= w830 <= 1
 0 <= w824 <= 1
 0 <= w823 <= 1
 0 <= w822 <= 1
 0 <= w821 <= 1
 0 <= w820 <= 1
 0 <= w814 <= 1
 0 <= w813 <= 1
 0 <= w812 <= 1
 0 <= w811 <= 1
 0 <= w810 <= 1
 0 <= w804 <= 1
 0 <= w803 <= 1
 0 <= w802 <= 1
 0 <= w801 <= 1
 0 <= w800 <= 1
 0 <= w744 <= 1
 0 <= w743 <= 1
 0 <= w742 <= 1
 0 <= w741 <= 1
 0 <= w740 <= 1
 0 <= w734 <= 1
 0 <= w733 <= 1
 0 <= w732 <= 1
 0 <= w731 <= 1
 0 <= w730 <= 1
 0 <= w724 <= 1
 0 <= w723 <= 1
 0 <= w722 <= 1
 0 <= w721 <= 1
 0 <= w720 <= 1
 0 <= w714 <= 1
 0 <= w713 <= 1
 0 <= w712 <= 1
 0 <= w711 <= 1
 0 <= w710 <= 1
 0 <= w704 <= 1
 0 <= w703 <= 1
 0 <= w702 <= 1
 0 <= w701 <= 1
 0 <= w700 <= 1
 0 <= w644 <= 1
 0 <= w643 <= 1
 0 <= w642 <= 1
 0 <= w641 <= 1
 0 <= w640 <= 1
 0 <= w634 <= 1
 0 <= w633 <= 1
 0 <= w632 <= 1
 0 <= w631 <= 1
 0 <= w630 <= 1
 0 <= w624 <= 1
 0 <= w623 <= 1
 0 <= w622 <= 1
 0 <= w621 <= 1
 0 <= w620 <= 1
 0 <= w614 <= 1
 0 <= w613 <= 1
 0 <= w612 <= 1
 0 <= w611 <= 1
 0 <= w610 <= 1
 0 <= w604 <= 1
 0 <= w603 <= 1
 0 <= w602 <= 1
 0 <= w601 <= 1
 0 <= w600 <= 1
 0 <= w544 <= 1
 0 <= w543 <= 1
 0 <= w542 <= 1
 0 <= w541 <= 1
 0 <= w540 <= 1
 0 <= w534 <= 1
 0 <= w533 <= 1
 0 <= w532 <= 1
 0 <= w531 <= 1
 0 <= w530 <= 1
 0 <= w524 <= 1
 0 <= w523 <= 1
 0 <= w522 <= 1
 0 <= w521 <= 1
 0 <= w520 <= 1
 0 <= w514 <= 1
 0 <= w513 <= 1
 0 <= w512 <= 1
 0 <= w511 <= 1
 0 <= w510 <= 1
 0 <= w504 <= 1
 0 <= w503 <= 1
 0 <= w502 <= 1
 0 <= w501 <= 1
 0 <= w500 <= 1
 0 <= w444 <= 1
 0 <= w443 <= 1
 0 <= w442 <= 1
 0 <= w441 <= 1
 0 <= w440 <= 1
 0 <= w434 <= 1
 0 <= w433 <= 1
 0 <= w432 <= 1
 0 <= w431 <= 1
 0 <= w430 <= 1
 0 <= w424 <= 1
 0 <= w423 <= 1
 0 <= w422 <= 1
 0 <= w421 <= 1
 0 <= w420 <= 1
 0 <= w414 <= 1
 0 <= w413 <= 1
 0 <= w412 <= 1
 0 <= w411 <= 1
 0 <= w410 <= 1
 0 <= w404 <= 1
 0 <= w403 <= 1
 0 <= w402 <= 1
 0 <= w401 <= 1
 0 <= w400 <= 1
 0 <= w344 <= 1
 0 <= w343 <= 1
 0 <= w342 <= 1
 0 <= w341 <= 1
 0 <= w340 <= 1
 0 <= w334 <= 1
 0 <= w333 <= 1
 0 <= w332 <= 1
 0 <= w331 <= 1
 0 <= w330 <= 1
 0 <= w324 <= 1
 0 <= w323 <= 1
 0 <= w322 <= 1
 0 <= w321 <= 1
 0 <= w320 <= 1
 0 <= w314 <= 1
 0 <= w313 <= 1
 0 <= w312 <= 1
 0 <= w311 <= 1
 0 <= w310 <= 1
 0 <= w304 <= 1
 0 <= w303 <= 1
 0 <= w302 <= 1
 0 <= w301 <= 1
 0 <= w300 <= 1
 0 <= w244 <= 1
 0 <= w243 <= 1
 0 <= w242 <= 1
 0 <= w241 <= 1
 0 <= w240 <= 1
 0 <= w234 <= 1
 0 <= w233 <= 1
 0 <= w232 <= 1
 0 <= w231 <= 1
 0 <= w230 <= 1
 0 <= w224 <= 1
 0 <= w223 <= 1
 0 <= w222 <= 1
 0 <= w221 <= 1
 0 <= w220 <= 1
 0 <= w214 <= 1
 0 <= w213 <= 1
 0 <= w212 <= 1
 0 <= w211 <= 1
 0 <= w210 <= 1
 0 <= w204 <= 1
 0 <= w203 <= 1
 0 <= w202 <= 1
 0 <= w201 <= 1
 0 <= w200 <= 1
 0 <= w144 <= 1
 0 <= w143 <= 1
 0 <= w142 <= 1
 0 <= w141 <= 1
 0 <= w140 <= 1
 0 <= w134 <= 1
 0 <= w133 <= 1
 0 <= w132 <= 1
 0 <= w131 <= 1
 0 <= w130 <= 1
 0 <= w124 <= 1
 0 <= w123 <= 1
 0 <= w122 <= 1
 0 <= w121 <= 1
 0 <= w120 <= 1
 0 <= w114 <= 1
 0 <= w113 <= 1
 0 <= w112 <= 1
 0 <= w111 <= 1
 0 <= w110 <= 1
 0 <= w104 <= 1
 0 <= w103 <= 1
 0 <= w102 <= 1
 0 <= w101 <= 1
 0 <= w100 <= 1
 0 <= w044 <= 1
 0 <= w043 <= 1
 0 <= w042 <= 1
 0 <= w041 <= 1
 0 <= w040 <= 1
 0 <= w034 <= 1
 0 <= w033 <= 1
 0 <= w032 <= 1
 0 <= w031 <= 1
 0 <= w030 <= 1
 0 <= w024 <= 1
 0 <= w023 <= 1
 0 <= w022 <= 1
 0 <= w021 <= 1
 0 <= w020 <= 1
 0 <= w014 <= 1
 0 <= w013 <= 1
 0 <= w012 <= 1
 0 <= w011 <= 1
 0 <= w010 <= 1
 0 <= w004 <= 1
 0 <= w003 <= 1
 0 <= w002 <= 1
 0 <= w001 <= 1
 0 <= w000 <= 1
 0 <= theta04 <= 1
 0 <= theta03 <= 1
 0 <= theta02 <= 1
 0 <= theta01 <= 1
 0 <= theta00 <= 1
 0 <= theta14 <= 1
 0 <= theta13 <= 1
 0 <= theta12 <= 1
 0 <= theta11 <= 1
 0 <= theta10 <= 1
 0 <= theta24 <= 1
 0 <= theta23 <= 1
 0 <= theta22 <= 1
 0 <= theta21 <= 1
 0 <= theta20 <= 1
 0 <= theta34 <= 1
 0 <= theta33 <= 1
 0 <= theta32 <= 1
 0 <= theta31 <= 1
 0 <= theta30 <= 1
 0 <= x00 <= 1
 0 <= y00 <= 1
 0 <= x01 <= 1
 0 <= y01 <= 1
 0 <= x02 <= 1
 0 <= y02 <= 1
 0 <= x03 <= 1
 0 <= y03 <= 1
 0 <= x04 <= 1
 0 <= y04 <= 1
 0 <= x10 <= 1
 0 <= y10 <= 1
 0 <= x11 <= 1
 0 <= y11 <= 1
 0 <= x12 <= 1
 0 <= y12 <= 1
 0 <= x13 <= 1
 0 <= y13 <= 1
 0 <= x14 <= 1
 0 <= y14 <= 1
 0 <= x20 <= 1
 0 <= y20 <= 1
 0 <= x21 <= 1
 0 <= y21 <= 1
 0 <= x22 <= 1
 0 <= y22 <= 1
 0 <= x23 <= 1
 0 <= y23 <= 1
 0 <= x24 <= 1
 0 <= y24 <= 1
 0 <= x30 <= 1
 0 <= y30 <= 1
 0 <= x31 <= 1
 0 <= y31 <= 1
 0 <= x32 <= 1
 0 <= y32 <= 1
 0 <= x33 <= 1
 0 <= y33 <= 1
 0 <= x34 <= 1
 0 <= y34 <= 1
 0 <= x40 <= 1
 0 <= y40 <= 1
 0 <= x41 <= 1
 0 <= y41 <= 1
 0 <= x42 <= 1
 0 <= y42 <= 1
 0 <= x43 <= 1
 0 <= y43 <= 1
 0 <= x44 <= 1
 0 <= y44 <= 1
 0 <= x50 <= 1
 0 <= y50 <= 1
 0 <= x51 <= 1
 0 <= y51 <= 1
 0 <= x52 <= 1
 0 <= y52 <= 1
 0 <= x53 <= 1
 0 <= y53 <= 1
 0 <= x54 <= 1
 0 <= y54 <= 1
 0 <= x60 <= 1
 0 <= y60 <= 1
 0 <= x61 <= 1
 0 <= y61 <= 1
 0 <= x62 <= 1
 0 <= y62 <= 1
 0 <= x63 <= 1
 0 <= y63 <= 1
 0 <= x64 <= 1
 0 <= y64 <= 1
 0 <= x70 <= 1
 0 <= y70 <= 1
 0 <= x71 <= 1
 0 <= y71 <= 1
 0 <= x72 <= 1
 0 <= y72 <= 1
 0 <= x73 <= 1
 0 <= y73 <= 1
 0 <= x74 <= 1
 0 <= y74 <= 1
 0 <= x80 <= 1
 0 <= y80 <= 1
 0 <= x81 <= 1
 0 <= y81 <= 1
 0 <= x82 <= 1
 0 <= y82 <= 1
 0 <= x83 <= 1
 0 <= y83 <= 1
 0 <= x84 <= 1
 0 <= y84 <= 1
 0 <= x90 <= 1
 0 <= y90 <= 1
 0 <= x91 <= 1
 0 <= y91 <= 1
 0 <= x92 <= 1
 0 <= y92 <= 1
 0 <= x93 <= 1
 0 <= y93 <= 1
 0 <= x94 <= 1
 0 <= y94 <= 1
 0 <= x100 <= 1
 0 <= y100 <= 1
 0 <= x101 <= 1
 0 <= y101 <= 1
 0 <= x102 <= 1
 0 <= y102 <= 1
 0 <= x103 <= 1
 0 <= y103 <= 1
 0 <= x104 <= 1
 0 <= y104 <= 1
 0 <= x110 <= 1
 0 <= y110 <= 1
 0 <= x111 <= 1
 0 <= y111 <= 1
 0 <= x112 <= 1
 0 <= y112 <= 1
 0 <= x113 <= 1
 0 <= y113 <= 1
 0 <= x114 <= 1
 0 <= y114 <= 1
 0 <= z0 <= 1
 0 <= z1 <= 1
 0 <= z2 <= 1
 0 <= z3 <= 1
 0 <= z4 <= 1
 0 <= z5 <= 1
 0 <= z6 <= 1
 0 <= z7 <= 1
 0 <= z8 <= 1
 0 <= z9 <= 1
 0 <= z10 <= 1
 0 <= z11 <= 1
End
