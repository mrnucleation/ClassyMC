  �  N   k820309    9          19.0        ��^                                                                                                          
       /home/SHARE/metastable/Classy/src/Script_TrajType.f90 INPUT_TRAJTYPE                                                     
                                                           
       GETXCOMMAND                                                     
       TRAJARRAY                   @               E               '�                    #TRAJ                �                                                        #TRAJECTORY                      @                              '                   #CLASSYCLASS    #DUMPFORCES %   #FILEUNIT &   #BOXNUM '   #OUTFREQ (   #FILENAME )   #SETUNIT *   #SETBOX .   #SETFILENAME 2   #SETFREQ 6   #OPENFILE :   #WRITEFRAME =   #CLOSEFILE A                �                                                          #CLASSYCLASS                      @                              '              
      #SCREENIO 	   #MAINTFREQ 
   #DESTRUCTOR    #EPILOGUE    #MAINTENANCE    #MODIFYIO    #PROLOGUE    #SAFETYCHECK    #SCREENOUT    #UPDATE "               �                               	                                          �                                                                         �                               
                                         �                                      
               10    1         �   �                       �                        #DESTRUCTOR    #         @                                                      #SELF              
                                                    #CLASSYCLASS    1         �   �                       �                        #EPILOGUE    #         @                                                      #SELF              
                                                    #CLASSYCLASS    1         �   �                       �                        #MAINTENANCE    #         @                                                      #SELF              
                                                    #CLASSYCLASS    1         �   �                       �                        #MODIFYIO    #         @                                                      #SELF    #LINE    #LINESTAT              
                                                    #CLASSYCLASS              
                                                    1                                                        1         �   �                       �                        #PROLOGUE    #         @                                                      #SELF              
                                                    #CLASSYCLASS    1         �   �                       �                        #SAFETYCHECK    #         @                                                      #SELF              
                                                    #CLASSYCLASS    1         �   �                       �                   	     #SCREENOUT     #         @                                                       #SELF !             
                                !                    #CLASSYCLASS    1         �   �                       �      "             
     #UPDATE #   #         @                                  #                    #SELF $             
                                $                    #CLASSYCLASS                �                               %                                         �                                                                         �                               &                                         �                                         ��������                        �                               '                                         �                                         ��������                        �                               (                                         �                                      �              5000                 �                              )     �                    1         �   �                       �      *              	    #SETUNIT +   #         @                                  +                    #SELF ,   #FILEUNIT -             
                                ,                   #TRAJECTORY              
                                  -           1         �   �                      �      .              
    #SETBOX /   #         @                                 /                    #SELF 0   #BOXNUM 1             
                                0                   #TRAJECTORY              
                                  1           1         �   �                      �      2             	     #SETFILENAME 3   #         @                                 3                    #SELF 4   #FILENAME 5             
                                4                   #TRAJECTORY              
                                5                    1 1         �   �                      �      6             
     #SETFREQ 7   #         @                                 7                    #SELF 8   #FREQ 9             
                                8                   #TRAJECTORY              
                                  9           1         �   �                      �      :                  #OPENFILE ;   #         @                                 ;                    #SELF <             
                                <                   #TRAJECTORY    1         �   �                       �      =                  #WRITEFRAME >   #         @                                  >                    #SELF ?   #ICYCLE @             
                                ?                   #TRAJECTORY              
                                 @           1         �   �                       �      A                  #CLOSEFILE B   #         @                                  B                    #SELF C             
                                 C                  #TRAJECTORY    #         @                                  D                    #LINE E   #COMMAND F   #COMNUM G   #LINESTAT H             
                                E                    1                                          F                     1           
                                  G                                                      H                     @ `                               I            �                        &                                           #TRJARRAY    #         @                                   J                    #LINE K   #TRAJNUM L   #LINESTAT M             
  @                             K                    1           
                                  L                     D @                               M               �   M      fn#fn    �   @   J   VARPRECISION    -  L   J  INPUT_FORMAT    y  J   J  TRAJDATA "   �  Z       TRJARRAY+TRAJDATA '     `   a   TRJARRAY%TRAJ+TRAJDATA .   }  
      TRAJECTORY+TRAJECTORYTEMPLATE :   �  a   a   TRAJECTORY%CLASSYCLASS+TRAJECTORYTEMPLATE +   �  �       CLASSYCLASS+MASTERTEMPLATE 4   �  �   a   CLASSYCLASS%SCREENIO+MASTERTEMPLATE 5   p  �   a   CLASSYCLASS%MAINTFREQ+MASTERTEMPLATE 6     X   a   CLASSYCLASS%DESTRUCTOR+MASTERTEMPLATE *   n  R       DESTRUCTOR+MASTERTEMPLATE /   �  Y   a   DESTRUCTOR%SELF+MASTERTEMPLATE 4     V   a   CLASSYCLASS%EPILOGUE+MASTERTEMPLATE (   o  R       EPILOGUE+MASTERTEMPLATE -   �  Y   a   EPILOGUE%SELF+MASTERTEMPLATE 7     Y   a   CLASSYCLASS%MAINTENANCE+MASTERTEMPLATE +   s  R       MAINTENANCE+MASTERTEMPLATE 0   �  Y   a   MAINTENANCE%SELF+MASTERTEMPLATE 4   	  V   a   CLASSYCLASS%MODIFYIO+MASTERTEMPLATE (   t	  j       MODIFYIO+MASTERTEMPLATE -   �	  Y   a   MODIFYIO%SELF+MASTERTEMPLATE -   7
  L   a   MODIFYIO%LINE+MASTERTEMPLATE 1   �
  @   a   MODIFYIO%LINESTAT+MASTERTEMPLATE 4   �
  V   a   CLASSYCLASS%PROLOGUE+MASTERTEMPLATE (     R       PROLOGUE+MASTERTEMPLATE -   k  Y   a   PROLOGUE%SELF+MASTERTEMPLATE 7   �  Y   a   CLASSYCLASS%SAFETYCHECK+MASTERTEMPLATE +     R       SAFETYCHECK+MASTERTEMPLATE 0   o  Y   a   SAFETYCHECK%SELF+MASTERTEMPLATE 5   �  W   a   CLASSYCLASS%SCREENOUT+MASTERTEMPLATE )     R       SCREENOUT+MASTERTEMPLATE .   q  Y   a   SCREENOUT%SELF+MASTERTEMPLATE 2   �  T   a   CLASSYCLASS%UPDATE+MASTERTEMPLATE &     R       UPDATE+MASTERTEMPLATE +   p  Y   a   UPDATE%SELF+MASTERTEMPLATE 9   �  �   a   TRAJECTORY%DUMPFORCES+TRAJECTORYTEMPLATE 7   m  �   a   TRAJECTORY%FILEUNIT+TRAJECTORYTEMPLATE 5     �   a   TRAJECTORY%BOXNUM+TRAJECTORYTEMPLATE 6   �  �   a   TRAJECTORY%OUTFREQ+TRAJECTORYTEMPLATE 7   ]  P   a   TRAJECTORY%FILENAME+TRAJECTORYTEMPLATE 6   �  U   a   TRAJECTORY%SETUNIT+TRAJECTORYTEMPLATE +     `       SETUNIT+TRAJECTORYTEMPLATE 0   b  X   a   SETUNIT%SELF+TRAJECTORYTEMPLATE 4   �  @   a   SETUNIT%FILEUNIT+TRAJECTORYTEMPLATE 5   �  T   a   TRAJECTORY%SETBOX+TRAJECTORYTEMPLATE *   N  ^       SETBOX+TRAJECTORYTEMPLATE /   �  X   a   SETBOX%SELF+TRAJECTORYTEMPLATE 1     @   a   SETBOX%BOXNUM+TRAJECTORYTEMPLATE :   D  Y   a   TRAJECTORY%SETFILENAME+TRAJECTORYTEMPLATE /   �  `       SETFILENAME+TRAJECTORYTEMPLATE 4   �  X   a   SETFILENAME%SELF+TRAJECTORYTEMPLATE 8   U  L   a   SETFILENAME%FILENAME+TRAJECTORYTEMPLATE 6   �  U   a   TRAJECTORY%SETFREQ+TRAJECTORYTEMPLATE +   �  \       SETFREQ+TRAJECTORYTEMPLATE 0   R  X   a   SETFREQ%SELF+TRAJECTORYTEMPLATE 0   �  @   a   SETFREQ%FREQ+TRAJECTORYTEMPLATE 7   �  V   a   TRAJECTORY%OPENFILE+TRAJECTORYTEMPLATE ,   @  R       OPENFILE+TRAJECTORYTEMPLATE 1   �  X   a   OPENFILE%SELF+TRAJECTORYTEMPLATE 9   �  X   a   TRAJECTORY%WRITEFRAME+TRAJECTORYTEMPLATE .   B  ^       WRITEFRAME+TRAJECTORYTEMPLATE 3   �  X   a   WRITEFRAME%SELF+TRAJECTORYTEMPLATE 5   �  @   a   WRITEFRAME%ICYCLE+TRAJECTORYTEMPLATE 8   8  W   a   TRAJECTORY%CLOSEFILE+TRAJECTORYTEMPLATE -   �  R       CLOSEFILE+TRAJECTORYTEMPLATE 2   �  X   a   CLOSEFILE%SELF+TRAJECTORYTEMPLATE )   9  y       GETXCOMMAND+INPUT_FORMAT .   �  L   a   GETXCOMMAND%LINE+INPUT_FORMAT 1   �  L   a   GETXCOMMAND%COMMAND+INPUT_FORMAT 0   J  @   a   GETXCOMMAND%COMNUM+INPUT_FORMAT 2   �  @   a   GETXCOMMAND%LINESTAT+INPUT_FORMAT #   �  �       TRAJARRAY+TRAJDATA     d  m       SCRIPT_TRAJTYPE %   �  L   a   SCRIPT_TRAJTYPE%LINE (     @   a   SCRIPT_TRAJTYPE%TRAJNUM )   ]  @   a   SCRIPT_TRAJTYPE%LINESTAT 