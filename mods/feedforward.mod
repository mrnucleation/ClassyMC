  �$  Q   k820309    l          18.0        ��]                                                                                                          
       ext/feedforward.f90 FEEDFORWARD              PRINT_MAT &         @                                      �                     #ARCH    #NETWORK              
                                                                  &                                           #         @                                                       #NET              
D                                      �              #NETWORK    #         @                                                       #NET    #FILE    #UNIT 	             
                                       �             #NETWORK              
 @                                                 1           
 @                               	           #         @                                   
                    #NET    #FILE    #UNIT              
                                       �             #NETWORK              
 @                                                 1           
 @                                          &         @                                      �                     #FILE    #UNIT    #NETWORK              
 @                                                 1           
 @                                          &         @                                      �                     #FILE    #UNIT    #NETWORK              
 @                                                 1           
 @                                          %         @                                                           #NET              
                                       �             #NETWORK    #         @                                                       #NET    #VERBOSE              
                                       �             #NETWORK              
 @                                          %         @                                                            #NET              
                                       �             #NETWORK    #         @                                                       #NET    #NW    #WEIGHTS              
D                                      �              #NETWORK              
                                                      
                                                      
    p          5 � p        r        5 � p        r                      #         @                                                       #NET               
D                                       �              #NETWORK    #         @                                   !                    #NET "   #NW #   #W_UPD $             
D                                 "     �              #NETWORK              
                                  #                    
                                  $                    
 	   p          5 � p        r #       5 � p        r #                     #         @                                   %                    #NET &   #NX '   #X (   #NY )   #VALUES *   #DERIVS ,   #Y -            
                                 &     �              #NETWORK              
                                  '                    
                                  (                    
 
   p          5 � p        r '       5 � p        r '                               
                                  )                    D @                               *                    
     p          5 8 �#NETWORK     p        r#NETWORK     &   U     +       5 8 �#NETWORK     p        r#NETWORK     &   U     +                              D @                               ,                    
     p          5 8 �#NETWORK     p        r#NETWORK     &   U     +       5 8 �#NETWORK     p        r#NETWORK     &   U     +                              D                                 -                    
     p          5 � p        r '       5 � p        r '                     #         @                                   .                    #NET /   #NX 0   #NY 1   #DERIVS 2   #JACOBIAN 3   #DY 5            
                                 /     �              #NETWORK              
                                  0                     
                                  1                    
                                  2                    
    p          5 8 �#NETWORK     p        r#NETWORK     /   U     +       5 8 �#NETWORK     p        r#NETWORK     /   U     +                              D                                 3                    
     p          5 8 �#NETWORK     p        r#NETWORK     /   U     4       5 8 �#NETWORK     p        r#NETWORK     /   U     4                              D                                 5                    
       p        5 � p        r 1   p          5 � p        r 1     5 � p        r 0       5 � p        r 1     5 � p        r 0                     #         @                                   6                    #NET 7   #NW 8   #NY 9   #VALUES :   #DERIVS ;   #JACOBIAN <   #DY_DW =            
                                 7     �              #NETWORK              
                                  8                     
                                  9                    
                                  :                    
    p          5 8 �#NETWORK     p        r#NETWORK     7   U     +       5 8 �#NETWORK     p        r#NETWORK     7   U     +                              
                                  ;                    
    p          5 8 �#NETWORK     p        r#NETWORK     7   U     +       5 8 �#NETWORK     p        r#NETWORK     7   U     +                              
                                  <                    
    p          5 8 �#NETWORK     p        r#NETWORK     7   U     4       5 8 �#NETWORK     p        r#NETWORK     7   U     4                              D                                 =                    
       p        5 � p        r 9   p          5 � p        r 9     5 � p        r 8       5 � p        r 9     5 � p        r 8                     #         @                                   >                    #NET ?   #ILAYER @   #FTYPE A             
D                                 ?     �              #NETWORK              
                                  @                     
  @                             A                    1 #         H                                 B                    #T C   #X D   #Y E   #DY F             
                                  C                     
  @                               D     
                D                                 E     
                 D                                 F     
                         @               A                '�                   #INIT G   #MEMSIZE H   #NLAYERS I   #NNODES J   #NNODES_MAX K   #F_A L   #W M   #IW N   #WSIZE 4   #IV O   #NVALUES +               �                               G                                          �                                                                          �                               H                               �                               I                             �                               J                                         &                                                       �                               K     X                        �                               L            `                             &                                                      �                               M            �                 
            &                                                      �                               N            �                             &                                                       �                               4     8      	                 �                               O            @             
               &                                                       �                               +     �               �   (      fn#fn !   �      b   uapp(FEEDFORWARD    �   g       NEW_NETWORK !   I  �   a   NEW_NETWORK%ARCH    �  Q       DEL_NETWORK     &  U   a   DEL_NETWORK%NET    {  e       SAVE_NETWORK !   �  U   a   SAVE_NETWORK%NET "   5  L   a   SAVE_NETWORK%FILE "   �  @   a   SAVE_NETWORK%UNIT #   �  e       SAVE_NETWORK_ASCII '   &  U   a   SAVE_NETWORK_ASCII%NET (   {  L   a   SAVE_NETWORK_ASCII%FILE (   �  @   a   SAVE_NETWORK_ASCII%UNIT      q       LOAD_NETWORK "   x  L   a   LOAD_NETWORK%FILE "   �  @   a   LOAD_NETWORK%UNIT #     q       LOAD_NETWORK_ASCII (   u  L   a   LOAD_NETWORK_ASCII%FILE (   �  @   a   LOAD_NETWORK_ASCII%UNIT      Y       FF_MEMSIZE    Z  U   a   FF_MEMSIZE%NET    �  ^       FF_PRINT_INFO "     U   a   FF_PRINT_INFO%NET &   b  @   a   FF_PRINT_INFO%VERBOSE     �  Y       FF_GET_NWEIGHTS $   �  U   a   FF_GET_NWEIGHTS%NET     P	  f       FF_INIT_WEIGHTS $   �	  U   a   FF_INIT_WEIGHTS%NET #   
  @   a   FF_INIT_WEIGHTS%NW (   K
  �   a   FF_INIT_WEIGHTS%WEIGHTS '   �
  Q       FF_RANDOM_INIT_WEIGHTS +   P  U   a   FF_RANDOM_INIT_WEIGHTS%NET "   �  d       FF_UPDATE_WEIGHTS &   	  U   a   FF_UPDATE_WEIGHTS%NET %   ^  @   a   FF_UPDATE_WEIGHTS%NW (   �  �   a   FF_UPDATE_WEIGHTS%W_UPD    R  �       FF_EVAL    �  U   a   FF_EVAL%NET    .  @   a   FF_EVAL%NX    n  �   a   FF_EVAL%X    "  @   a   FF_EVAL%NY    b    a   FF_EVAL%VALUES    n    a   FF_EVAL%DERIVS    z  �   a   FF_EVAL%Y    .  �       FF_DERIV    �  U   a   FF_DERIV%NET      @   a   FF_DERIV%NX    F  @   a   FF_DERIV%NY     �    a   FF_DERIV%DERIVS "   �    a   FF_DERIV%JACOBIAN    �  $  a   FF_DERIV%DY    �  �       FF_WDERIV    T  U   a   FF_WDERIV%NET    �  @   a   FF_WDERIV%NW    �  @   a   FF_WDERIV%NY !   )    a   FF_WDERIV%VALUES !   5    a   FF_WDERIV%DERIVS #   A    a   FF_WDERIV%JACOBIAN     M  $  a   FF_WDERIV%DY_DW %   q  h       FF_CHANGE_ACTIVATION )   �  U   a   FF_CHANGE_ACTIVATION%NET ,   .  @   a   FF_CHANGE_ACTIVATION%ILAYER +   n  L   a   FF_CHANGE_ACTIVATION%FTYPE    �  e       FF_ACTIVATE      @   a   FF_ACTIVATE%T    _  @   a   FF_ACTIVATE%X    �  @   a   FF_ACTIVATE%Y    �  @   a   FF_ACTIVATE%DY      �       NETWORK    �  �   a   NETWORK%INIT     �   H   a   NETWORK%MEMSIZE     �   H   a   NETWORK%NLAYERS    !  �   a   NETWORK%NNODES #   �!  H   a   NETWORK%NNODES_MAX    �!  �   a   NETWORK%F_A    �"  �   a   NETWORK%W    #  �   a   NETWORK%IW    �#  H   a   NETWORK%WSIZE    �#  �   a   NETWORK%IV     �$  H   a   NETWORK%NVALUES 