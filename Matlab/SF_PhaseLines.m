function out=Phaselines(material_1, material_2, varargin)
    %v. 1.0 Clinton & Journaux, 2020
    %Function calculates phase equilibria coordinates (K,MPa) based on the 
    %Gibbs energy of two materials. 
    %
    %Option 'Plot' develops a figure showcasing calculated data and 'meta'
    %adds predicted metstable extensions.
    %
    %Example:
    % >> Phaselines('Ih', 'water1', 'plot', 'meta')
    
    %Creation of varargin options
    options = [0 0];
    if (~isempty(varargin))
        for c=1:length(varargin)
            switch  varargin{c}
                case {'meta'}
                    options(1)=1;
                case {'plot'}
                    options(2)=1;
                otherwise
                    error('Invalid argument. Check format and available options.')
            end
        end
    end
    
    if  options(1)==1
        %Pressure and Temperature Structures for Extension
        
        PhaselinesD.Ih.water1_Ih = {(0.1:300),(239:0.5:285)};
        PhaselinesD.Ih.II_Ih = {(0.1:240), (0.1:0.5:269)};
        PhaselinesD.Ih.III_Ih = {(200:215), (230:0.5:260)};
        PhaselinesD.II.III_II = {(100:495), (218.5:0.5:258.2)};
        PhaselinesD.II.V_II = {(204:807), (164:0.5:270)};
        PhaselinesD.II.VI_II = {(0.1:900), (0.1:0.5:270)};
        PhaselinesD.III.V_III = {(340:391),(210:0.5:270)};
        PhaselinesD.III.water1_III = {(47:490), (240:0.5:260)};
        PhaselinesD.V.water1_V ={(179:725),(239:0.5:290)};
        PhaselinesD.V.VI_V ={(633:684),(122:0.5:297)};
        PhaselinesD.VI.water1_VI ={(337:2300),(270:0.5:375)};
        
    else
        %Pressure and Temperature without Metastable Extensions
        PhaselinesD.Ih.water1_Ih = {(0.7:207.7),(250:0.5:273)};
        PhaselinesD.Ih.II_Ih = {(1:239), (1:0.5:239)};
        PhaselinesD.Ih.III_Ih = {(207.5:210.5), (235:0.5:255)};
        PhaselinesD.II.III_II = {(205:356), (235:0.5:251)};
        PhaselinesD.II.V_II = {(354:671), (190:0.5:250)};
        PhaselinesD.II.VI_II = {(1:900), (1:0.5:202)};
        PhaselinesD.III.V_III = {(350:356),(249.3:0.5:257.3)};
        PhaselinesD.III.water1_III = {(205:0.5:350.5), (250:0.5:257)};
        PhaselinesD.V.water1_V ={(340:635),(255:0.5:275)};
        PhaselinesD.V.VI_V ={(630:680),(201.9:0.5:274)};
        PhaselinesD.VI.water1_VI ={(630:2300),(273.38:0.5:375)};
    end
    
    
    %Option: Meta and Plotting
    if options(1)==1 && options(2)==1
        ylabel('Temperature [K]');
        xlabel('Pressure [MPa]');
        hold on
        %Data 'Ih' and 'water1'
        if ((strcmp('Ih', material_1) && strcmp('water1', material_2)) || strcmp('water1', material_1) && strcmp('Ih', material_2))
            A=SeaFreeze(PhaselinesD.Ih.water1_Ih, 'Ih');
            B=SeaFreeze(PhaselinesD.Ih.water1_Ih, 'water1');
            Gibbs1= A.G;
            Gibbs2= B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.Ih.water1_Ih{1}(1,:);
            T=PhaselinesD.Ih.water1_Ih{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            hold on
            plot(207.5930,251.1191,'bo')
            hold on
            TP(:,1:40)=[];
            TP(:,64)=[251.1191,207.5930];
            plot(TP(2,1:64), TP(1,1:64),':r');
            hold on
            plot(TP(2,64:316), TP(1,64:316), '-r');
            %Data for 'Ih' and 'II'
        elseif ((strcmp('Ih', material_1) && strcmp('II', material_2)) || strcmp('II', material_1) && strcmp('Ih', material_2))
            A = SeaFreeze(PhaselinesD.Ih.II_Ih, 'Ih');
            B = SeaFreeze(PhaselinesD.Ih.II_Ih, 'II');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.Ih.II_Ih{1}(1,:);
            T=PhaselinesD.Ih.II_Ih{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,1)=[];
            TP(:,722:773)=[];
            TP(:,685)=[238.2371,209.885];
            plot(209.885,238.2371,'bo')
            hold on
            plot(TP(2,236:685), TP(1,236:685), '-r');
            hold on
            plot(TP(2,685:721), TP(1, 685:721), ':r')
            hold on
            plot(TP(2,1:236), TP(1,1:236), 'LineStyle', '--', 'Color', 'r')
            %Data 'Ih' and 'III'
        elseif((strcmp('Ih', material_1) && strcmp('III', material_2)) || strcmp('III', material_1) && strcmp('Ih', material_2))
            A = SeaFreeze(PhaselinesD.Ih.III_Ih, 'Ih');
            B = SeaFreeze(PhaselinesD.Ih.III_Ih, 'III');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.Ih.III_Ih{1}(1,:);
            T=PhaselinesD.Ih.III_Ih{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,1)=[];
            TP(:,20)=[238.2377, 209.885];
            TP(:,47)=[251.1191,207.5930];
            plot(209.885,238.2377, 'bo')
            hold on
            plot(207.5930,251.1191, 'bo')
            hold on
            plot(TP(2,20:47), TP(1,20:47), '-r');
            hold on
            plot(TP(2,1:20), TP(1, 1:20), ':r');
            plot(TP(2,47:66), TP(1, 47:66), ':r');
            %Data for 'II' and 'III'
        elseif ((strcmp('II', material_1) && strcmp('III', material_2)) || strcmp('III', material_1) && strcmp('II', material_2))
            A = SeaFreeze(PhaselinesD.II.III_II, 'II');
            B = SeaFreeze(PhaselinesD.II.III_II, 'III');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.II.III_II{1}(1,:);
            T=PhaselinesD.II.III_II{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,1)=[];
            TP(:,323)=[238.2379,209.885];
            TP(:,155)=[249.4176,355.5042];
            plot(355.5042,249.4176, 'bo')
            hold on
            plot(209.885,238.2379,'bo')
            hold on
            plot(TP(2,155:323), TP(1,155:323), '-r');
            hold on
            plot(TP(2,1:155),TP(1,1:155), ':r');
            hold on
            plot(TP(2,323:453),TP(1,323:453), ':r');
            %Data for 'II' and 'V'
        elseif ((strcmp('II', material_1) && strcmp('V', material_2)) || strcmp('V', material_1) && strcmp('II', material_2))
            A = SeaFreeze(PhaselinesD.II.V_II, 'II');
            B = SeaFreeze(PhaselinesD.II.V_II, 'V');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.II.V_II{1}(1,:);
            T=PhaselinesD.II.V_II{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,1)=[];
            TP(:,592)=[249.4176,355.5042];
            TP(:,181)=[201.9335,670.8401];
            plot(670.8401,201.9335, 'bo')
            hold on
            plot(355.5042,249.4176, 'bo')
            hold on
            plot(TP(2,181:592), TP(1,181:592), '-r');
            hold on
            plot(TP(2,1:181), TP(1,1:181), ':r');
            hold on
            plot(TP(2,592:643), TP(1,592:643), ':r');
            TP(:,643:784)= [];
            %Data for 'II' and 'VI'
        elseif ((strcmp('II', material_1) && strcmp('VI', material_2)) || strcmp('VI', material_1) && strcmp('II', material_2))
            A = SeaFreeze(PhaselinesD.II.VI_II, 'II');
            B = SeaFreeze(PhaselinesD.II.VI_II, 'VI');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.II.VI_II{1}(1,:);
            T=PhaselinesD.II.VI_II{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,1)=[];
            TP(:,459)=[201.9335,670.8401];
            plot(670.8401,201.9335, 'bo')
            hold on
            plot(TP(2,459:786),TP(1,459:786),':r');
            TP(:,786:835)=[];
            hold on
            % II/VI Add on
            load('iceII_sp_G_fPT_1200MPa.mat', 'sp_G_fPT');
            P = 898:1200;
            T = 0.1:0.5:90;
            PhaselinesD.II.VI_II_add={P,T};
            out=SeaFreeze(PhaselinesD.II.VI_II_add, 'VI');
            materialI=out.G;
            out=fnGval(sp_G_fPT, PhaselinesD.II.VI_II_add);
            materialII=out.G;
            Z= materialI-materialII;
            A=contourc(T,P,Z, [0,0]);
            A(:,1:3)=[];
            A(:,1)=[86.7641,899.100];
            %Combining A and PT for add on
            A=fliplr(A);
            TP=cat(2,A,TP);
            hold on
            plot(TP(2,1:505),TP(1,1:505), '--r');
            hold on
            plot(TP(2,505:806),TP(1,505:806),'-r');
            %Data for 'III' and 'V'
        elseif ((strcmp('III', material_1) && strcmp('V', material_2)) || strcmp('V', material_1) && strcmp('III', material_2))
            A = SeaFreeze(PhaselinesD.III.V_III, 'III');
            B = SeaFreeze(PhaselinesD.III.V_III, 'V');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.III.V_III{1}(1,:);
            T=PhaselinesD.III.V_III{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,1)=[];
            TP(:,133)=[256.1641,350.1095];
            TP(:,115)=[249.4176,355.5042];
            plot(355.5042,249.4176,'bo');
            hold on
            plot(350.1095,256.1641,'bo');
            hold on
            plot(TP(2,115:133), TP(1,115:133), '-r');
            hold on
            plot(TP(2,1:115), TP(1,1:115),':r');
            hold on
            plot(TP(2,133:171), TP(1,133:171), ':r')
            %Data for 'III' and 'water1'
        elseif ((strcmp('III', material_1) && strcmp('water1', material_2)) || strcmp('III', material_1) && strcmp('water1', material_2))
            A = SeaFreeze(PhaselinesD.III.water1_III, 'III');
            B = SeaFreeze(PhaselinesD.III.water1_III, 'water1');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.III.water1_III{1}(1,:);
            T=PhaselinesD.III.water1_III{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,1)=[];
            TP(:,146)=[256.1641,350.1095];
            TP(:,299)=[251.1191, 207.5930];
            plot(350.1095,256.1641, 'bo')
            hold on
            plot(207.5930,251.1191, 'bo')
            hold on
            plot(TP(2,146:299), TP(1,146:299), '-r');
            hold on
            plot(TP(2,1:146), TP(1,1:146), ':r');
            hold on
            plot(TP(2,299:477),TP(1,299:477), ':r');
            %Data for 'V' and 'water1'
        elseif ((strcmp('V', material_1) && strcmp('water1', material_2)) || strcmp('water1', material_1) && strcmp('V', material_2))
            A = SeaFreeze(PhaselinesD.V.water1_V, 'V');
            B = SeaFreeze(PhaselinesD.V.water1_V, 'water1');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.V.water1_V{1}(1,:);
            T=PhaselinesD.V.water1_V{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,1)=[];
            TP(:,100)=[273.40653,634.3997];
            TP(:,418)=[256.1641,350.1095];
            plot(634.3997,273.40653, 'bo')
            hold on
            plot(350.1095,256.1641, 'bo')
            hold on
            plot(TP(2,100:418), TP(1,100:418), '-r');
            hold on
            plot(TP(2,1:100), TP(1,1:100), ':r');
            hold on
            plot(TP(2,418:575), TP(1,418:575), ':r');
            TP(:,576:622)=[];
            %Data for 'VI' and 'water1'
        elseif ((strcmp('VI', material_1) && strcmp('water1', material_2)) || strcmp('water1', material_1) && strcmp('VI', material_2))
            A = SeaFreeze(PhaselinesD.VI.water1_VI, 'VI');
            B = SeaFreeze(PhaselinesD.VI.water1_VI, 'water1');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.VI.water1_VI{1}(1,:);
            T=PhaselinesD.VI.water1_VI{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,1)=[];
            TP(:,1835)=[273.4063,634.399];
            plot(634.399,273.4063, 'bo')
            hold on
            plot(TP(2,1:1835), TP(1,1:1835), '-r');
            hold on
            plot(TP(2,1821:1879), TP(1,1821:1879), ':r');
            %Data for 'V' and 'VI'
        elseif ((strcmp('V', material_1) && strcmp('VI', material_2)) || strcmp('VI', material_1) && strcmp('V', material_2))
            A = SeaFreeze(PhaselinesD.V.VI_V, 'V');
            B = SeaFreeze(PhaselinesD.V.VI_V, 'VI');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.V.VI_V{1}(1,:);
            T=PhaselinesD.V.VI_V{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,1)=[];
            TP(:,302)=[273.4063,634.39];
            TP(:,111)=[201.9335,670.8409];
            plot(634.39,273.4063,'bo')
            hold on
            plot(670.8409,201.9335,'bo')
            hold on
            plot(TP(2,111:302), TP(1,111:302), '-r');
            hold on
            plot(TP(2,1:111), TP(1,1:111), ':r');
            hold on
            plot(TP(2,302:305), TP(1,302:305), ':r');
        else
            disp('Error: Please type materials and options as ''(''x'',''y'',''plot'',''meta'')'' or ''(''x'',''y'',''meta'',''plot'')')
        end
        
        % Option: Plotting without Metastable Extensions
    elseif options(2)==1 && options(1)==0
        ylabel('Temperature [K]');
        xlabel('Pressure [MPa]');
        hold on
        %Data for 'Ih' and 'water1'
        if ((strcmp('Ih', material_1) && strcmp('water1', material_2)) || strcmp('water1', material_1) && strcmp('Ih', material_2))
            A=SeaFreeze(PhaselinesD.Ih.water1_Ih, 'Ih');
            B=SeaFreeze(PhaselinesD.Ih.water1_Ih, 'water1');
            Gibbs1= A.G;
            Gibbs2= B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.Ih.water1_Ih{1}(1,:);
            T=PhaselinesD.Ih.water1_Ih{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,4)=[251.1191;207.5930];
            TP(:,1:3)=[];
            plot(TP(2,:),TP(1,:), '-r')
            %Data for 'Ih' and 'II'
        elseif ((strcmp('Ih', material_1) && strcmp('II', material_2)) || strcmp('II', material_1) && strcmp('Ih', material_2))
            A = SeaFreeze(PhaselinesD.Ih.II_Ih, 'Ih');
            B = SeaFreeze(PhaselinesD.Ih.II_Ih, 'II');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.Ih.II_Ih{1}(1,:);
            T=PhaselinesD.Ih.II_Ih{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,657)=[238.2371,209.885];
            TP(:,658:659)=[];
            TP(:,1)=[];
            plot(TP(2,206:656), TP(1,206:656), '-r');
            hold on
            plot(TP(2,1:206), TP(1,1:206), 'Linestyle', '--', 'Color', 'r')
            %Data for 'Ih' and 'III'
        elseif((strcmp('Ih', material_1) && strcmp('III', material_2)) || strcmp('III', material_1) && strcmp('Ih', material_2))
            A = SeaFreeze(PhaselinesD.Ih.III_Ih, 'Ih');
            B = SeaFreeze(PhaselinesD.Ih.III_Ih, 'III');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.Ih.III_Ih{1}(1,:);
            T=PhaselinesD.Ih.III_Ih{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,8)=[238.2377,209.885];
            TP(:,37)=[251.1191, 207.5930];
            TP(:,38)=[];
            TP(:,1:7)=[];
            plot(TP(2,:), TP(1,:), '-r');
            %Data for 'II' and 'III'
        elseif ((strcmp('II', material_1) && strcmp('III', material_2)) || strcmp('III', material_1) && strcmp('II', material_2))
            A = SeaFreeze(PhaselinesD.II.III_II, 'II');
            B = SeaFreeze(PhaselinesD.II.III_II, 'III');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.II.III_II{1}(1,:);
            T=PhaselinesD.II.III_II{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,2)=[249.4176, 355.5042];
            TP(:,169)=[238.2379,209.885];
            TP(:,170:176)= [];
            TP(:,1)=[];
            plot(TP(2,:), TP(1,:), '-r');
            %Data for 'II' and 'V'
        elseif ((strcmp('II', material_1) && strcmp('V', material_2)) || strcmp('V', material_1) && strcmp('II', material_2))
            A = SeaFreeze(PhaselinesD.II.V_II, 'II');
            B = SeaFreeze(PhaselinesD.II.V_II, 'V');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.II.V_II{1}(1,:);
            T=PhaselinesD.II.V_II{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,1)=[];
            TP(:,413:414)=[];
            TP(:,1)=[201.9335,670.8401];
            TP(:,412)=[249.4176,355.5042];
            plot(TP(2,:), TP(1,:), '-r');
            %Data for 'II' and 'VI'
        elseif ((strcmp('II', material_1) && strcmp('VI', material_2)) || strcmp('VI', material_1) && strcmp('II', material_2))
            A = SeaFreeze(PhaselinesD.II.VI_II, 'II');
            B = SeaFreeze(PhaselinesD.II.VI_II, 'VI');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.II.VI_II{1}(1,:);
            T=PhaselinesD.II.VI_II{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,463)=[201.9335,670.8401];
            TP(:,1)=[];
            hold on
            %Add on to VI and II
            load('iceII_sp_G_fPT_1200MPa.mat', 'sp_G_fPT');
            P = 899:1200;
            T = 0.1:0.5:90;
            PhaselinesD.II.VI_II_add={P,T};
            out=SeaFreeze(PhaselinesD.II.VI_II_add, 'VI');
            materialI=out.G;
            out=fnGval(sp_G_fPT, PhaselinesD.II.VI_II_add);
            materialII=out.G;
            Z= materialI-materialII;
            A=contourc(T,P,Z, [0,0]);
            A(:,1)=[];
            A(:,1)=[86.7641,899.100];
            %Combining the Addition
            A=fliplr(A);
            TP = cat(2,A,TP);
            hold on
            plot(TP(2,1:505),TP(1,1:505),'--r');
            hold on
            plot(TP(2,505:809),TP(1,505:809),'-r');
            %Data for 'III' and 'V'
        elseif ((strcmp('III', material_1) && strcmp('V', material_2)) || strcmp('V', material_1) && strcmp('III', material_2))
            A = SeaFreeze(PhaselinesD.III.V_III, 'III');
            B = SeaFreeze(PhaselinesD.III.V_III, 'V');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.III.V_III{1}(1,:);
            T=PhaselinesD.III.V_III{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,2) = [249.4176,355.5042];
            TP(:,20)=[256.1641,350.1095];
            TP(:,21:22)=[];
            TP(:,1)= [];
            plot(TP(2,:), TP(1,:), '-r');
            %Data for 'III' and 'water1'
        elseif ((strcmp('III', material_1) && strcmp('water1', material_2)) || strcmp('water1', material_1) && strcmp('III', material_2))
            A = SeaFreeze(PhaselinesD.III.water1_III, 'III');
            B = SeaFreeze(PhaselinesD.III.water1_III, 'water1');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.III.water1_III{1}(1,:);
            T=PhaselinesD.III.water1_III{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,3)=[256.1641,350.1095];
            TP(:,297)=[251.1191,207.5930];
            TP(:,298:304)=[];
            TP(:,1:2)=[];
            plot(TP(2,:), TP(1,:), '-r');
            %Data for 'V' and 'water1'
        elseif ((strcmp('V', material_1) && strcmp('water1', material_2)) || strcmp('water1', material_1) && strcmp('V', material_2))
            A = SeaFreeze(PhaselinesD.V.water1_V, 'V');
            B = SeaFreeze(PhaselinesD.V.water1_V, 'water1');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.V.water1_V{1}(1,:);
            T=PhaselinesD.V.water1_V{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,2)=[273.4063,634.3997];
            TP(:,319)=[256.1641,350.1095];
            TP(:,320:333)=[];
            TP(:,1)=[];
            plot(TP(2,:), TP(1,:), '-r');
            %Data for 'VI' and 'water1'
        elseif ((strcmp('VI', material_1) && strcmp('water1', material_2)) || strcmp('water1', material_1) && strcmp('VI', material_2))
            A = SeaFreeze(PhaselinesD.VI.water1_VI, 'VI');
            B = SeaFreeze(PhaselinesD.VI.water1_VI, 'water1');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.VI.water1_VI{1}(1,:);
            T=PhaselinesD.VI.water1_VI{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,1834)=[273.4063,634.3997];
            TP(:,1835)=[];
            TP(:,1)=[];
            plot(TP(2,:), TP(1,:), '-r');
            %Data for 'V' and 'VI'
        elseif ((strcmp('V', material_1) && strcmp('VI', material_2)) || strcmp('VI', material_1) && strcmp('V', material_2))
            A = SeaFreeze(PhaselinesD.V.VI_V, 'V');
            B = SeaFreeze(PhaselinesD.V.VI_V, 'VI');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.V.VI_V{1}(1,:);
            T=PhaselinesD.V.VI_V{2}(1,:);
            [TP,~] = contour(T, P, Z, [0, 0], 'Visible', 'off');
            TP(:,194)=[273.4063,634.3997];
            TP(:,2)=[201.9335,670.8401];
            TP(:,195)=[];
            TP(:,1)=[];
            plot(TP(2,:), TP(1,:), '-r');
        else
            disp('Error: Please type materials and option as ''(''x'', ''y'', ''plot'')')
        end
        
        %Option: Without Plotting or Metastable Ext.
    elseif options(1)== 0  && options(2)==0
        %Data for 'Ih' and 'water1'
        if ((strcmp('Ih', material_1) && strcmp('water1', material_2)) || strcmp('water1', material_1) && strcmp('Ih', material_2))
            A=SeaFreeze(PhaselinesD.Ih.water1_Ih, 'Ih');
            B=SeaFreeze(PhaselinesD.Ih.water1_Ih, 'water1');
            Gibbs1= A.G;
            Gibbs2= B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.Ih.water1_Ih{1}(1,:);
            T=PhaselinesD.Ih.water1_Ih{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,2)=[251.1191;207.5930];
            TP(:,1)=[];
            %Data for 'Ih' and 'II'
        elseif ((strcmp('Ih', material_1) && strcmp('II', material_2)) || strcmp('II', material_1) && strcmp('Ih', material_2))
            A = SeaFreeze(PhaselinesD.Ih.II_Ih, 'Ih');
            B = SeaFreeze(PhaselinesD.Ih.II_Ih, 'II');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.Ih.II_Ih{1}(1,:);
            T=PhaselinesD.Ih.II_Ih{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,656)=[238.2377,209.885];
            TP(:,657:659)=[];
            TP(:,1)=[];
            %Data for 'Ih' and 'III'
        elseif((strcmp('Ih', material_1) && strcmp('III', material_2)) || strcmp('III', material_1) && strcmp('Ih', material_2))
            A = SeaFreeze(PhaselinesD.Ih.III_Ih, 'Ih');
            B = SeaFreeze(PhaselinesD.Ih.III_Ih, 'III');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.Ih.III_Ih{1}(1,:);
            T=PhaselinesD.Ih.III_Ih{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,32)=[238.2377,209.885];
            TP(:,13)=[251.1191, 207.5930];
            TP(:,33:38)=[];
            TP(:,1:12)=[];
            %Data for 'II' and 'III'
        elseif ((strcmp('II', material_1) && strcmp('III', material_2)) || strcmp('III', material_1) && strcmp('II', material_2))
            A = SeaFreeze(PhaselinesD.II.III_II, 'II');
            B = SeaFreeze(PhaselinesD.II.III_II, 'III');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.II.III_II{1}(1,:);
            T=PhaselinesD.II.III_II{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,176)=[249.4176, 355.5042];
            TP(:,7)=[238.2377,209.885];
            TP(:,1:6)=[];
            %Data for 'II' and 'V'
        elseif ((strcmp('II', material_1) && strcmp('V', material_2)) || strcmp('V', material_1) && strcmp('II', material_2))
            A = SeaFreeze(PhaselinesD.II.V_II, 'II');
            B = SeaFreeze(PhaselinesD.II.V_II, 'V');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.II.V_II{1}(1,:);
            T=PhaselinesD.II.V_II{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,1)=[];
            TP(:,413:414)=[];
            TP(:,1)=[201.9335,670.8401];
            TP(:,412)=[249.4176,355.5042];
            %Data for 'II' and 'VI'
        elseif ((strcmp('II', material_1) && strcmp('VI', material_2)) || strcmp('VI', material_1) && strcmp('II', material_2))
            A = SeaFreeze(PhaselinesD.II.VI_II, 'II');
            B = SeaFreeze(PhaselinesD.II.VI_II, 'VI');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.II.VI_II{1}(1,:);
            T=PhaselinesD.II.VI_II{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,463)=[201.9335,670.8401];
            TP(:,1)=[];
            %Data for 'III' and 'V'
        elseif ((strcmp('III', material_1) && strcmp('V', material_2)) || strcmp('V', material_1) && strcmp('III', material_2))
            A = SeaFreeze(PhaselinesD.III.V_III, 'III');
            B = SeaFreeze(PhaselinesD.III.V_III, 'V');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.III.V_III{1}(1,:);
            T=PhaselinesD.III.V_III{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,22) = [249.4176,355.5042];
            TP(:,3)=[256.1641,350.1095];
            TP(:,1:2)= [];
            %Data for 'III' and 'water1'
        elseif ((strcmp('III', material_1) && strcmp('water1', material_2)) || strcmp('water1', material_1) && strcmp('III', material_2))
            A = SeaFreeze(PhaselinesD.III.water1_III, 'III');
            B = SeaFreeze(PhaselinesD.III.water1_III, 'water1');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.III.water1_III{1}(1,:);
            T=PhaselinesD.III.water1_III{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,300)=[256.1641,350.1095];
            TP(:,8)=[251.1191,207.5930];
            TP(:,301:304)=[];
            TP(:,1:7)=[];
            %Data for 'V' and 'water1'
        elseif ((strcmp('V', material_1) && strcmp('water1', material_2)) || strcmp('water1', material_1) && strcmp('V', material_2))
            A = SeaFreeze(PhaselinesD.V.water1_V, 'V');
            B = SeaFreeze(PhaselinesD.V.water1_V, 'water1');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.V.water1_V{1}(1,:);
            T=PhaselinesD.V.water1_V{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,333)=[273.4063,634.3997];
            TP(:,14)=[256.1641,350.1095];
            TP(:,1:13)=[];
            %Data for 'VI' and 'water1'
        elseif ((strcmp('VI', material_1) && strcmp('water1', material_2)) || strcmp('water1', material_1) && strcmp('VI', material_2))
            A = SeaFreeze(PhaselinesD.VI.water1_VI, 'VI');
            B = SeaFreeze(PhaselinesD.VI.water1_VI, 'water1');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.VI.water1_VI{1}(1,:);
            T=PhaselinesD.VI.water1_VI{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,1834)=[273.4063,634.3997];
            TP(:,1835)=[];
            TP(:,1)=[];
            %Data for 'V' and 'VI'
        elseif ((strcmp('V', material_1) && strcmp('VI', material_2)) || strcmp('VI', material_1) && strcmp('V', material_2))
            A = SeaFreeze(PhaselinesD.V.VI_V, 'V');
            B = SeaFreeze(PhaselinesD.V.VI_V, 'VI');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.V.VI_V{1}(1,:);
            T=PhaselinesD.V.VI_V{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,194)=[273.4063,634.3997];
            TP(:,2)=[201.9335,670.8401];
            TP(:,195)=[];
            TP(:,1)=[];
        else
            disp('Error: Please type materials as ''x', 'y''')
        end
        
        %Option: Meta only
    elseif options(1)==1
        %Data for 'Ih' and 'water1'
        if((strcmp('Ih', material_1) && strcmp('water1', material_2)) || strcmp('water1', material_1) && strcmp('Ih', material_2))
            A=SeaFreeze(PhaselinesD.Ih.water1_Ih, 'Ih');
            B=SeaFreeze(PhaselinesD.Ih.water1_Ih, 'water1');
            Gibbs1= A.G;
            Gibbs2= B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.Ih.water1_Ih{1}(1,:);
            T=PhaselinesD.Ih.water1_Ih{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,1)=[];
            TP(:,317:355)=[];
            %Data for 'Ih' and 'II'
        elseif ((strcmp('Ih', material_1) && strcmp('II', material_2)) || strcmp('II', material_1) && strcmp('Ih', material_2))
            A = SeaFreeze(PhaselinesD.Ih.II_Ih, 'Ih');
            B = SeaFreeze(PhaselinesD.Ih.II_Ih, 'II');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.Ih.II_Ih{1}(1,:);
            T=PhaselinesD.Ih.II_Ih{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,1)=[];
            TP(:,722:773)=[];
            %Data for 'Ih' and 'III'
        elseif((strcmp('Ih', material_1) && strcmp('III', material_2)) || strcmp('III', material_1) && strcmp('Ih', material_2))
            A = SeaFreeze(PhaselinesD.Ih.III_Ih, 'Ih');
            B = SeaFreeze(PhaselinesD.Ih.III_Ih, 'III');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.Ih.III_Ih{1}(1,:);
            T=PhaselinesD.Ih.III_Ih{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,1)=[];
            %Data for 'II' and 'III'
        elseif ((strcmp('II', material_1) && strcmp('III', material_2)) || strcmp('III', material_1) && strcmp('II', material_2))
            A = SeaFreeze(PhaselinesD.II.III_II, 'II');
            B = SeaFreeze(PhaselinesD.II.III_II, 'III');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.II.III_II{1}(1,:);
            T=PhaselinesD.II.III_II{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,1)=[];
            %Data for 'II' and 'V'
        elseif ((strcmp('II', material_1) && strcmp('V', material_2)) || strcmp('V', material_1) && strcmp('II', material_2))
            A = SeaFreeze(PhaselinesD.II.V_II, 'II');
            B = SeaFreeze(PhaselinesD.II.V_II, 'V');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.II.V_II{1}(1,:);
            T=PhaselinesD.II.V_II{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,1)=[];
            TP(:,643:784)=[];
            %Data for 'II' and 'VI'
        elseif ((strcmp('II', material_1) && strcmp('VI', material_2)) || strcmp('VI', material_1) && strcmp('II', material_2))
            A = SeaFreeze(PhaselinesD.II.VI_II, 'II');
            B = SeaFreeze(PhaselinesD.II.VI_II, 'VI');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.II.VI_II{1}(1,:);
            T=PhaselinesD.II.VI_II{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,1)=[];
            TP(:,786:835)=[];
                        %Data for 'III' and 'V'
        elseif ((strcmp('III', material_1) && strcmp('V', material_2)) || strcmp('V', material_1) && strcmp('III', material_2))
            A = SeaFreeze(PhaselinesD.III.V_III, 'III');
            B = SeaFreeze(PhaselinesD.III.V_III, 'V');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.III.V_III{1}(1,:);
            T=PhaselinesD.III.V_III{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,1)=[];
            %Data for 'III' and 'water1'
        elseif ((strcmp('III', material_1) && strcmp('water1', material_2)) || strcmp('water1', material_1) && strcmp('III', material_2))
            A = SeaFreeze(PhaselinesD.III.water1_III, 'III');
            B = SeaFreeze(PhaselinesD.III.water1_III, 'water1');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.III.water1_III{1}(1,:);
            T=PhaselinesD.III.water1_III{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,1)=[];
            %Data for 'V' and 'water1'
        elseif ((strcmp('V', material_1) && strcmp('water1', material_2)) || strcmp('water1', material_1) && strcmp('V', material_2))
            A = SeaFreeze(PhaselinesD.V.water1_V, 'V');
            B = SeaFreeze(PhaselinesD.V.water1_V, 'water1');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.V.water1_V{1}(1,:);
            T=PhaselinesD.V.water1_V{2}(1,:);
            TP= contourc(T, P, Z, [0, 0]);
            TP(:,1)=[];
            TP(:,1:47)=[];
            %Data for 'VI' and 'water1'
        elseif ((strcmp('VI', material_1) && strcmp('water1', material_2)) || strcmp('water1', material_1) && strcmp('VI', material_2))
            A = SeaFreeze(PhaselinesD.VI.water1_VI, 'VI');
            B = SeaFreeze(PhaselinesD.VI.water1_VI, 'water1');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.VI.water1_VI{1}(1,:);
            T=PhaselinesD.VI.water1_VI{2}(1,:);
            TP= contourc(T, P, Z, [0, 0]);
            TP(:,1)=[];
            %Data for 'V' and 'VI'
        elseif ((strcmp('V', material_1) && strcmp('VI', material_2)) || strcmp('VI', material_1) && strcmp('V', material_2))
            A = SeaFreeze(PhaselinesD.V.VI_V, 'V');
            B = SeaFreeze(PhaselinesD.V.VI_V, 'VI');
            Gibbs1=A.G;
            Gibbs2=B.G;
            Z=Gibbs1-Gibbs2;
            P=PhaselinesD.V.VI_V{1}(1,:);
            T=PhaselinesD.V.VI_V{2}(1,:);
            TP = contourc(T, P, Z, [0, 0]);
            TP(:,1)=[];
        else
            disp('Error: Please type materials and option as ''(''x'',''y'',''meta'')')
        end
    end
    out=flip(TP',2);
    
