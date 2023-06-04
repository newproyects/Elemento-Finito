% draw a finite element model in the current figure
%

cla; % Clear Axes
figNumber=gcf; % Crea una ventana emergente para mostrar la figura
set(figNumber,'renderer','painters'); % Fija propiedades como 'renderer', 'painters' a la figura
txtHndl=get(figNumber,'UserData'); % Consulta propiedades antes fijadas, información

tic

num_nodes = size(coords, 1); % 6 es el número de nodos
num_elems = size(elem_con, 2); % 2 regiones cuadradas encerradas por los 6 nodos 
num_bc_vals = size(bc_vals, 1); % 6 es el tamaño de los vectores bc_vals
num_bcs = size(bcs, 1); % 6 es el tamaño de vectores bcs
xmax = 10; xmin = 0; ymax = 10; ymin = 0; zmax = 0.1; zmin = 0; % Rango para los ejes
ax=[xmin xmax ymin ymax]; % Variable que fija el rango para ejes en 2D

% Se ajusta para un tamaño optimizado de la imagen en función de los extremos de la figura.
% Después se busca el valor más pequeño dentro de los nodos disponibles en
% sus coordenadas para reajustar

if(num_nodes > 0)
    xmin = min(coords(:, 1)); 
    xmax = max(coords(:, 1));
    ymin = min(coords(:, 2));
    ymax = max(coords(:, 2));
    zmin = min(coords(:, 3));
    zmax = max(coords(:, 3));
    ax_size = [max( 2, (xmax - xmin)) max( 2, (ymax - ymin)) (zmax - zmin)]; 
    mag_ax_size = norm(ax_size);
    xmin = xmin - 0.2*ax_size(1);
    xmax = xmax + 0.2*ax_size(1);
    ymin = ymin - 0.2*ax_size(2);
    ymax = ymax + 0.2*ax_size(2);
    zmin = zmin - 0.2*ax_size(3);
    zmax = zmax + 0.2*ax_size(3);
    if(abs(zmin - zmax) > .0000001) %% Si existe una dimensión Z entonces...
        ax=[xmin xmax ymin ymax zmin zmax];
    else                            %% Caso contrario dibujemos en 2D
        ax=[xmin xmax ymin ymax];        
    end
    
    % draw bc_vals
    % 
    
    if(draw_bc_vals == 1)
        trans_bc_max_mag = max( max( abs( bc_vals( 1:num_bc_vals, 1:3 ) ) ) ); % Escalar Valor máximo de las fuerzas de traslación 
        rot_bc_max_mag = max(max(abs(bc_vals(1:num_bc_vals, 4:6))));  % Escalar Valor máximo de las fuerzas de rotación (torcas)
        
        %%
        
        for(i = 1 : num_bc_vals); %% Para todos los vectores en bc_vals
            for(dir = 1 : 3) %% Para los componentes de dichos vectores
                if(bc_vals(i, dir) ~= 0) %% Cuando dichos componentes sean diferentes de cero
                    load_vec = zeros(3, 1); %% Creamos un vector receptor en 3D con entradas nulas
                    load_vec(dir) = bc_vals(i, dir)/trans_bc_max_mag; %% Normalizamos los componentes del vector
                    
                    line([coords(i, 1) coords(i, 1) + load_vec(1)*mag_ax_size/10],... 
                        [coords(i, 2) coords(i, 2) + load_vec(2)*mag_ax_size/10],...
                        [coords(i, 3) coords(i, 3) + load_vec(3)*mag_ax_size/10],...
                        'LineStyle','-', ...
                        'Marker','none', ...
                        'Color','black', ...
                        'MarkerSize',25, ...
                        'HitTest', 'off'); %% 
                    
                    loc = [coords(i, 1) + load_vec(1)*mag_ax_size/10,...
                           coords(i, 2) + load_vec(2)*mag_ax_size/10,...
                           coords(i, 3) + load_vec(3)*mag_ax_size/10];
                    text(loc(1), loc(2), loc(3), num2str(abs(bc_vals(i, dir))), 'color', 'g');
                end
            end
            for(dir = 1 : 3)
                if(bc_vals(i, dir + 3) ~= 0)
                    load_vec = zeros(3, 1);
                    load_vec(dir) = bc_vals(i, dir + 3)/rot_bc_max_mag;
                    line([coords(i, 1) coords(i, 1) + load_vec(1)*mag_ax_size/7],...
                        [coords(i, 2) coords(i, 2) + load_vec(2)*mag_ax_size/7],...
                        [coords(i, 3) coords(i, 3) + load_vec(3)*mag_ax_size/7],...
                        'LineStyle','-', ...
                        'Marker','none', ...
                        'Color','r', ...
                        'MarkerSize',25, ...
                        'HitTest', 'off');
                    loc = [coords(i, 1) + load_vec(1)*mag_ax_size/7,...
                           coords(i, 2) + load_vec(2)*mag_ax_size/7,...
                           coords(i, 3) + load_vec(3)*mag_ax_size/7];
                    text(loc(1), loc(2), loc(3), num2str(bc_vals(i, dir+3)), 'color', 'r');
                end
            end
        end
    end
    % draw bcs
    if(draw_bcs == 1)
        for(i = 1 : num_bcs);
            for(dir = 1 :3)
                if(bcs(i, dir) ~= 0)
                    bc_vec = zeros(3, 1);
                    bc_vec(dir) = bcs(i, dir);
                    line([coords(i, 1) coords(i, 1) + bc_vec(1)*mag_ax_size/30],...
                         [coords(i, 2) coords(i, 2) + bc_vec(2)*mag_ax_size/30],...
                         [coords(i, 3) coords(i, 3) + bc_vec(3)*mag_ax_size/30],...
                        'LineStyle','-', ...
                        'LineWidth', 4, ...
                        'Marker','none', ...
                        'color',[.3 .3 .3], ...
                        'MarkerSize',25, ...                        
                        'HitTest', 'off');
                end
            end

            for(dir = 1 : 3)
                if(bcs(i, dir + 3) ~= 0)
                    bc_vec = zeros(3, 1);
                    bc_vec(dir) = bcs(i, dir + 3);
                    line([coords(i, 1) coords(i, 1) + bc_vec(1)*mag_ax_size/30],...
                         [coords(i, 2) coords(i, 2) + bc_vec(2)*mag_ax_size/30],...
                         [coords(i, 3) coords(i, 3) + bc_vec(3)*mag_ax_size/30],...
                        'LineStyle','-', ...
                        'LineWidth', 2, ...
                        'Marker','none', ...
                        'color','r', ...
                        'MarkerSize',25, ...
                        'HitTest', 'off');
                end
            end

        end
    end
    %
    
    % draw nodes
    tic
    if(draw_nodes == 1)
        for(i = 1 : num_nodes);
            line(coords(i, 1), coords(i, 2), coords(i, 3), ...
                'LineStyle', 'none', ...
                'Marker', '.', ...
                'Color', 'b', ...
                'MarkerSize', 20, ...                
                'UserData', [1 i],...
                'ButtonDownFcn','perifem(''hitobject'')');
            if(draw_node_nums == 1)
                text(coords(i, 1), coords(i, 2), coords(i, 3), [' ' num2str(i)], 'color', 'b',...
                    'UserData', [1 i], 'ButtonDownFcn','perifem(''hitobject'')');
            end
        end
    end
    
end






    
% draw elements

if (norm(draw_elems) ~= 0)
    for i = 1 : num_elems
            if (mode(1:1) == 'a')
                draw_elem(i, scale_factor, shrink_factor, 'hit_off', 'black',elem_props, nod_disps, draw_elems, coords, elem_con, ...
                    draw_elem_nums, view_max_st, draw_def_elems);
            else
                draw_elem(i, scale_factor, shrink_factor, 'hit_onn', 'black', elem_props,nod_disps, draw_elems, coords, elem_con, ...
                    draw_elem_nums, view_max_st, draw_def_elems);
            end
    end
end

axis(ax);
view(2);
hold on;
set(gca,...
    'Box','on',...
    'UserData',[],...
    'DataAspectRatio', [1 1 1],...
    'ButtonDownFcn','perifem(''ptselect'')');
xlabel('X');
ylabel('Y');
zlabel('Z');

toc


return