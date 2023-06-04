function draw_elem(elem_num, scale_factor, shrink_fact, hit, color, elem_props, nod_disps, draw_elems, coords, elem_con, ...
    draw_elem_nums, view_max_st, draw_def_elems)
%draw elements, both deformed and undeformed
%


elem_type = elem_props{elem_num}(1);
%
% Draw deformed elements first
%
if(~isempty(nod_disps))


    if(elem_type == 1 & (draw_def_elems(1) == 1))
        d1 = nod_disps(elem_con{elem_num}(1), :)*scale_factor;
        d2 = nod_disps(elem_con{elem_num}(2), :)*scale_factor;
        d1 = coords(elem_con{elem_num}(1), :) + d1(1:3);
        d2 = coords(elem_con{elem_num}(2), :) + d2(1:3);
        linewidth = 2*elem_type;
        line([d1(1) d2(1)],...
            [d1(2) d2(2)],...
            [d1(3) d2(3)],...
            'LineStyle','-', ...
            'LineWidth', linewidth, ...
            'Marker','none', ...
            'Color', [0 0 0], ...
            'MarkerSize',25, ...
            'EraseMode','none',...
            'HitTest', 'off');
    elseif(elem_type == 6 & (draw_def_elems(6) == 1))
        linewidth = 9;
       
        xo = coords(elem_con{elem_num}(1), 1) + nod_disps(elem_con{elem_num}(1), 1)*scale_factor;
        yo = coords(elem_con{elem_num}(1), 2) + nod_disps(elem_con{elem_num}(1), 2)*scale_factor;
        x = xo + (20*shrink_fact)*cos([0:10:360]*(pi/180));
        y = yo + (20*shrink_fact)*sin([0:10:360]*(pi/180));

        if(hit(1:7) == 'hit_off')
            line(x,y,...
                'LineStyle','--', ...
                'LineWidth', 0.1*linewidth, ...
                'Marker','none', ...
                'Color', color, ...
                'MarkerSize',25, ...
                'EraseMode','none',...
                'HitTest', 'off');
        else
            line(x,y,...
                'LineStyle','-', ...
                'LineWidth', 0.1*linewidth, ...
                'Marker','none', ...
                'Color', color, ...
                'MarkerSize',25, ...
                'EraseMode','none',...
                'HitTest', 'off');
        end

    elseif(elem_type == 5 & (draw_def_elems(5) == 1))

        linewidth = 2;
        c(1, :) = (coords(elem_con{elem_num}(1), :) + scale_factor*nod_disps(elem_con{elem_num}(1), 1:3))*shrink_fact + (coords(elem_con{elem_num}(2), :) + scale_factor*nod_disps(elem_con{elem_num}(2), 1:3))*(1 - shrink_fact);
        c(2, :) = (coords(elem_con{elem_num}(1), :) + scale_factor*nod_disps(elem_con{elem_num}(1), 1:3))*(1 - shrink_fact) + (coords(elem_con{elem_num}(2), :) + scale_factor*nod_disps(elem_con{elem_num}(2), 1:3))*shrink_fact;
        if(hit(1:7) == 'hit_off')
            line(c(:, 1),...
                c(:, 2),...
                c(:, 3),...
                'LineStyle','-', ...
                'LineWidth', linewidth, ...
                'Marker','none', ...
                'Color', [0.1 0.1 0.1], ...
                'MarkerSize',25, ...
                'EraseMode','none',...
                'HitTest', 'off');
        else
            line(c(:, 1),...
                c(:, 2),...
                c(:, 3),...
                'LineStyle','-', ...
                'LineWidth', linewidth, ...
                'Marker','none', ...
                'Color', [0.1 0.1 0.1], ...
                'MarkerSize',25, ...
                'EraseMode','none',...
                'UserData', [2 elem_num],...
                'ButtonDownFcn','perifem(''hitobject'')');
        end


    elseif(elem_type == 2 & (draw_def_elems(2) == 1))
        L = norm(coords(elem_con{elem_num}(1), :) - coords(elem_con{elem_num}(2), :));
        disp_vec = [nod_disps(elem_con{elem_num}(1), :)  nod_disps(elem_con{elem_num}(2), :)]';
        elem_coords = [coords(elem_con{elem_num}(1), :); coords(elem_con{elem_num}(2), :)];
        Rotation_Angle = elem_props{elem_num}(8);
        T = frame_trans_matrix(Rotation_Angle, elem_coords);
        local_disp_vec = T*disp_vec;
        num_segs = 10;
        for i = 1 : num_segs
            xi(1:2) = [(i-1)/num_segs i/num_segs];
            NL1(1:2) = 1-xi;
            NL2(1:2) = xi;
            N1(1:2) = 2*xi.^3 - 3*xi.^2 + 1;
            N2(1:2) = L*(xi.^3 - 2*xi.^2 + xi);
            N3(1:2) =  -2*xi.^3 + 3*xi.^2;
            N4(1:2) = L*(xi.^3 - xi.^2);
            shape_1 = [NL1(1) 0     0     0  0     0      NL2(1) 0     0     0  0     0;...
                0       N1(1) 0     0  0     N2(1)  0      N3(1) 0     0  0     N4(1);...
                0       0     N1(1) 0  -N2(1) 0      0      0     N3(1) 0  -N4(1) 0];

            shape_2 = [NL1(2) 0     0     0  0     0      NL2(2) 0     0     0   0     0;...
                0       N1(2) 0     0  0     N2(2)  0      N3(2) 0     0  0     N4(2);...
                0       0     N1(2) 0  -N2(2) 0      0      0    N3(2) 0  -N4(2) 0];
            local_disp1 = shape_1*local_disp_vec;
            local_disp2 = shape_2*local_disp_vec;
            disp1 = T(1:3, 1:3)'*local_disp1;
            disp2 = T(1:3, 1:3)'*local_disp2;
            d1 = NL1(1)*elem_coords(1, :) + NL2(1)*elem_coords(2, :)  + disp1(1:3)'*scale_factor;
            d2 = NL1(2)*elem_coords(1, :) + NL2(2)*elem_coords(2, :)  + disp2(1:3)'*scale_factor;
            linewidth = 2*elem_type;
            line([d1(1) d2(1)],...
                [d1(2) d2(2)],...
                [d1(3) d2(3)],...
                'LineStyle','-', ...
                'LineWidth', linewidth, ...
                'Marker','none', ...
                'Color', [.5 .5 .5], ...
                'MarkerSize',25, ...
                'EraseMode','none',...
                'HitTest', 'off');
        end

    elseif(elem_type == 3 & (draw_def_elems(3) == 1))
        elem_num_nodes = elem_props{elem_num}(2);
        linewidth = 1;
        for node = 1 : elem_num_nodes
            el_coords(node, :) = coords(elem_con{elem_num}(node), :);
            el_disps(node, :) = nod_disps(elem_con{elem_num}(node), 1:3);
        end
        for(node = 1 : elem_num_nodes)
            if(elem_num_nodes == 3)
                if(node == 1)
                    xi = [1, 0];
                elseif(node == 2)
                    xi = [0, 1];
                elseif(node == 3)
                    xi = [0, 0];
                end
            elseif(elem_num_nodes == 4)
                if(node == 1)
                    xi = [-1, -1];
                elseif(node == 2)
                    xi = [1, -1];
                elseif(node == 3)
                    xi = [1, 1];
                elseif(node == 4)
                    xi = [-1, 1];
                end
            elseif(elem_num_nodes == 6)
                if(node == 1)
                    xi = [0, 0];
                elseif(node == 2)
                    xi = [0.5, 0];
                elseif(node == 3)
                    xi = [1, 0];
                elseif(node == 4)
                    xi = [0.5, 0.5];
                elseif(node == 5)
                    xi = [0, 1];
                elseif(node == 6)
                    xi = [0, 0.5];
                end
            elseif(elem_num_nodes == 8)
                if(node == 1)
                    xi = [-1, -1];
                elseif(node == 2)
                    xi = [0, -1];
                elseif(node == 3)
                    xi = [1, -1];
                elseif(node == 4)
                    xi = [1, 0];
                elseif(node == 5)
                    xi = [1, 1];
                elseif(node == 6)
                    xi = [0, 1];
                elseif(node == 7)
                    xi = [-1, 1];
                elseif(node == 8)
                    xi = [-1, 0];
                end
            end
            N = shape_functions(xi, 'membrane', elem_num_nodes);
            c(node, :) = N*(el_coords + scale_factor*el_disps);
        end
        c(elem_num_nodes + 1, :) = c(1, :);
        line(c(:, 1),c(:, 2),c(:, 3),...
            'LineStyle','-', ...
            'LineWidth', linewidth, ...
            'Marker','none', ...
            'Color', 'black', ...
            'MarkerSize',25, ...            
            'HitTest', 'off');
    elseif(elem_type == 4 & (draw_def_elems(3) == 1))
        elem_num_nodes = elem_props{elem_num}(2);
        linewidth = 1;
        for node = 1 : elem_num_nodes
            el_coords(node, :) = coords(elem_con{elem_num}(node), :);
            el_disps(node, :) = nod_disps(elem_con{elem_num}(node), 1:3);
        end
        for(node = 1 : elem_num_nodes)
            if(elem_num_nodes == 4)
                if(node == 1)
                    xi = [1, 0, 0];
                elseif(node == 2)
                    xi = [0, 1, 0];
                elseif(node == 3)
                    xi = [0, 0, 1];
                elseif(node == 4)
                    xi = [0, 0, 0];
                end
            elseif(elem_num_nodes == 10)
                if(node == 1)
                    xi = [1, 0, 0];
                elseif(node == 2)
                    xi = [.5, .5, 0];
                elseif(node == 3)
                    xi = [0, 1, 0];
                elseif(node == 4)
                    xi = [0, .5, .5];
                elseif(node == 5)
                    xi = [0, 0, 1];
                elseif(node == 6)
                    xi = [.5, 0, .5];
                elseif(node == 7)
                    xi = [.5, 0, 0];
                elseif(node == 8)
                    xi = [0, .5, 0];
                elseif(node == 9)
                    xi = [0, 0, .5];
                elseif(node == 10)
                    xi = [0, 0, 0];
                end
            elseif(elem_num_nodes == 8)
                if(node == 1)
                    xi = [-1, -1, -1];
                elseif(node == 2)
                    xi = [1, -1, -1];
                elseif(node == 3)
                    xi = [1, 1, -1];
                elseif(node == 4)
                    xi = [-1, 1, -1];
                elseif(node == 5)
                    xi = [-1, -1, 1];
                elseif(node == 6)
                    xi = [1, -1, 1];
                elseif(node == 7)
                    xi = [1, 1, 1];
                elseif(node == 8)
                    xi = [-1, 1, 1];
                end
            elseif(elem_num_nodes == 20)
                if(node == 1)
                    xi = [-1, -1, -1];
                elseif(node == 2)
                    xi = [0, -1, -1];
                elseif(node == 3)
                    xi = [1, -1, -1];
                elseif(node == 4)
                    xi = [1, 0, -1];
                elseif(node == 5)
                    xi = [1, 1, -1];
                elseif(node == 6)
                    xi = [0, 1, -1];
                elseif(node == 7)
                    xi = [-1, 1, -1];
                elseif(node == 8)
                    xi = [-1, 0, -1];
                elseif(node == 9)
                    xi = [-1, -1, 0];
                elseif(node == 10)
                    xi = [1, -1, 0];
                elseif(node == 11)
                    xi = [1, 1, 0];
                elseif(node == 12)
                    xi = [-1, 1, 0];
                elseif(node == 13)
                    xi = [-1, -1, 1];
                elseif(node == 14)
                    xi = [0, -1, 1];
                elseif(node == 15)
                    xi = [1, -1, 1];
                elseif(node == 16)
                    xi = [1, 0, 1];
                elseif(node == 17)
                    xi = [1, 1, 1];
                elseif(node == 18)
                    xi = [0, 1, 1];
                elseif(node == 19)
                    xi = [-1, 1, -1];
                elseif(node == 20)
                    xi = [-1, 0, -1];
                end
            end
            N = shape_functions(xi, 'solid', elem_num_nodes);
            c(node, :) = N*(el_coords + scale_factor*el_disps);
        end

        if(elem_num_nodes == 4)
            edges = [c(1,:); c(2,:); c(3,:); c(1,:); c(4,:); c(2,:); c(3,:); c(4,:)];
            xi = [0.25 0.25 0.25];
        elseif(elem_num_nodes == 10)
            xi = [0.25 0.25 0.25];
        elseif(elem_num_nodes == 8)
            edges = [c(1,:); c(2,:); c(3,:); c(4,:); c(1,:);...
                c(5,:); c(6,:); c(7,:); c(8,:); c(5,:);...
                c(6,:); c(2,:); c(3,:); c(7,:); c(8,:); c(4,:)];
            xi = [0 0 0];
        elseif(elem_num_nodes == 20)
            xi = [0 0 0];
        end

        line(edges(:, 1), edges(:, 2), edges(:, 3),...
            'LineStyle','-', ...
            'LineWidth', linewidth, ...
            'Marker','none', ...
            'Color',[.5 .5 .5], ...
            'MarkerSize',25, ...
             'HitTest', 'off');
    else

    end
end
%
% Draw undeformed elements
%
if((elem_type <= 2) & (draw_elems(elem_type)==1))
    linewidth = 2*elem_type;
    c(1, :) = coords(elem_con{elem_num}(1), :)*shrink_fact + coords(elem_con{elem_num}(2), :)*(1 - shrink_fact);
    c(2, :) = coords(elem_con{elem_num}(1), :)*(1 - shrink_fact) + coords(elem_con{elem_num}(2), :)*shrink_fact;
    if(hit(1:7) == 'hit_off')
        line(c(:, 1),...
            c(:, 2),...
            c(:, 3),...
            'LineStyle','-', ...
            'LineWidth', linewidth, ...
            'Marker','none', ...
            'Color', color, ...
            'MarkerSize',25, ...
            'EraseMode','none',...
            'HitTest', 'off');
    else
        line(c(:, 1),...
            c(:, 2),...
            c(:, 3),...
            'LineStyle','-', ...
            'LineWidth', linewidth, ...
            'Marker','none', ...
            'Color', color, ...
            'MarkerSize',25, ...
            'EraseMode','none',...
            'UserData', [2 elem_num],...
            'ButtonDownFcn','perifem(''hitobject'')');
    end
    if(draw_elem_nums == 1)
        mid = ( coords(elem_con{elem_num}(1), :) + coords(elem_con{elem_num}(2), :) )/2;
        text(mid(1), mid(2), mid(3), [' ' num2str(elem_num)], 'color', 'b',...
            'UserData', [2 elem_num], 'ButtonDownFcn', 'perifem(''hitobject'')');
    end
    if(hit(1:7) == 'hit_off')
        line(c(:, 1),c(:, 2),c(:, 3),...
            'LineStyle','-', ...
            'LineWidth', linewidth, ...
            'Marker','none', ...
            'Color', color, ...
            'MarkerSize',25, ...
            'EraseMode','none',...
            'HitTest', 'off');
    else
        line(c(:, 1),c(:, 2),c(:, 3),...
            'LineStyle','-', ...
            'LineWidth', linewidth, ...
            'Marker','none', ...
            'Color', color, ...
            'MarkerSize',25, ...
            'EraseMode','none',...
            'UserData', [2 elem_num],...
            'ButtonDownFcn','perifem(''hitobject'')');
    end
    if(draw_elem_nums == 1)
        text(mid(1), mid(2), mid(3), [' ' num2str(elem_num)], 'color', 'g',...
            'UserData', [2 elem_num], 'ButtonDownFcn','perifem(''hitobject'')');
    end

elseif ((elem_type == 5)& (draw_elems(elem_type)==1))

    linewidth = 2;
    c(1, :) = coords(elem_con{elem_num}(1), :)*shrink_fact + coords(elem_con{elem_num}(2), :)*(1 - shrink_fact);
    c(2, :) = coords(elem_con{elem_num}(1), :)*(1 - shrink_fact) + coords(elem_con{elem_num}(2), :)*shrink_fact;
    if(hit(1:7) == 'hit_off')
        line(c(:, 1),...
            c(:, 2),...
            c(:, 3),...
            'LineStyle','-', ...
            'LineWidth', linewidth, ...
            'Marker','none', ...
            'Color', [0.1 0.1 0.1], ...
            'MarkerSize',25, ...
            'EraseMode','none',...
            'HitTest', 'off');
    else
        line(c(:, 1),...
            c(:, 2),...
            c(:, 3),...
            'LineStyle','-', ...
            'LineWidth', linewidth, ...
            'Marker','none', ...
            'Color', [0.1 0.1 0.1], ...
            'MarkerSize',25, ...
            'EraseMode','none',...
            'UserData', [2 elem_num],...
            'ButtonDownFcn','perifem(''hitobject'')');
    end
    if(draw_elem_nums == 1)
        mid = ( coords(elem_con{elem_num}(1), :) + coords(elem_con{elem_num}(2), :) )/2;
        text(mid(1), mid(2), mid(3), [' ' num2str(elem_num)], 'color', 'b',...
            'UserData', [2 elem_num], 'ButtonDownFcn', 'perifem(''hitobject'')');
    end
    if(hit(1:7) == 'hit_off')
        line(c(:, 1),c(:, 2),c(:, 3),...
            'LineStyle','-', ...
            'LineWidth', linewidth, ...
            'Marker','none', ...
            'Color', [0.2 0.2 0.2], ...
            'MarkerSize',25, ...
            'EraseMode','none',...
            'HitTest', 'off');
    else
        line(c(:, 1),c(:, 2),c(:, 3),...
            'LineStyle','-', ...
            'LineWidth', linewidth, ...
            'Marker','none', ...
            'Color', [0.7 0.7 0.7], ...
            'MarkerSize',25, ...
            'EraseMode','none',...
            'UserData', [2 elem_num],...
            'ButtonDownFcn','perifem(''hitobject'')');
    end
    if(draw_elem_nums == 1)
        text(mid(1), mid(2), mid(3), [' ' num2str(elem_num)], 'color', 'g',...
            'UserData', [2 elem_num], 'ButtonDownFcn','perifem(''hitobject'')');
    end

elseif ((elem_type == 6)& (draw_elems(elem_type)==1))


    if(hit(1:7) == 'hit_off')
        linewidth = 2;
        xo = coords(elem_con{elem_num}(1), 1);
        yo = coords(elem_con{elem_num}(1), 2);
        x = xo + (10*shrink_fact)*cos([0:10:360]*(pi/180));
        y = yo + (10*shrink_fact)*sin([0:10:360]*(pi/180));
        line(x,y,...
            'LineStyle','-', ...
            'LineWidth', linewidth, ...
            'Marker','none', ...
            'Color', [0.8 0.8 0.8], ...
            'MarkerSize',5, ...
            'EraseMode','none',...
            'HitTest', 'off');
    else
        linewidth = 2;
        xo = coords(elem_con{elem_num}(1), 1);
        yo = coords(elem_con{elem_num}(1), 2);
        x = xo + (10*shrink_fact)*cos([0:10:360]*(pi/180));
        y = yo + (10*shrink_fact)*sin([0:10:360]*(pi/180));
        line(x,y,...
            'LineStyle','-', ...
            'LineWidth', linewidth, ...
            'Marker','none', ...
            'Color', [0.8 0.8 0.8], ...
            'MarkerSize',25, ...
            'EraseMode','none',...
            'UserData', [2 elem_num],...
            'ButtonDownFcn','perifem(''hitobject'')');
    end
    if(draw_elem_nums == 1)
        text(xo + (shrink_fact), yo, [' ' num2str(elem_num)], 'color', 'g',...
            'UserData', [2 elem_num], 'ButtonDownFcn','perifem(''hitobject'')');
    end

elseif(elem_type == 3 & (draw_elems(elem_type)==1))
    elem_num_nodes = elem_props{elem_num}(2);
    linewidth = 1;
    for node = 1 : elem_num_nodes
        el_coords(node, :) = coords(elem_con{elem_num}(node), :);
    end
    for(node = 1 : elem_num_nodes)
        if(elem_num_nodes == 3)
            if(node == 1)
                xi = [1-shrink_fact, shrink_fact/2];
            elseif(node == 2)
                xi = [shrink_fact/2, 1-shrink_fact];
            elseif(node == 3)
                xi = [shrink_fact/2, shrink_fact/2];
            end
        elseif(elem_num_nodes == 4)
            if(node == 1)
                xi = [-1+shrink_fact, -1+shrink_fact];
            elseif(node == 2)
                xi = [1-shrink_fact, -1+shrink_fact];
            elseif(node == 3)
                xi = [1-shrink_fact, 1-shrink_fact];
            elseif(node == 4)
                xi = [-1+shrink_fact, 1-shrink_fact];
            end
        elseif(elem_num_nodes == 6)
            if(node == 1)
                xi = [shrink_fact, shrink_fact];
            elseif(node == 2)
                xi = [0.5, shrink_fact];
            elseif(node == 3)
                xi = [1-2.4142*shrink_fact, shrink_fact];
            elseif(node == 4)
                xi = [0.5 - 0.7071*shrink_fact, 0.5 - 0.7071*shrink_fact];
            elseif(node == 5)
                xi = [shrink_fact, 1 - 2.4142*shrink_fact];
            elseif(node == 6)
                xi = [shrink_fact, 0.5];
            end
        elseif(elem_num_nodes == 8)
            if(node == 1)
                xi = [-1+shrink_fact, -1+shrink_fact];
            elseif(node == 2)
                xi = [0, -1+shrink_fact];
            elseif(node == 3)
                xi = [1-shrink_fact, -1+shrink_fact];
            elseif(node == 4)
                xi = [1-shrink_fact, 0];
            elseif(node == 5)
                xi = [1-shrink_fact, 1-shrink_fact];
            elseif(node == 6)
                xi = [0, 1-shrink_fact];
            elseif(node == 7)
                xi = [-1+shrink_fact, 1-shrink_fact];
            elseif(node == 8)
                xi = [-1+shrink_fact, 0];
            end
        end
        N = shape_functions(xi, 'membrane', elem_num_nodes);
        c(node, :) = N*el_coords;
    end
    c(elem_num_nodes + 1, :) = c(1, :);
    if(elem_num_nodes == 3)
        xi = [0.33333 0.33333];
    elseif(elem_num_nodes == 4)
        xi = [0 0];
    elseif(elem_num_nodes == 6)
        xi = [0.33333 0.33333];
    elseif(elem_num_nodes == 8)
        xi = [0 0];
    end
    N = shape_functions(xi, 'membrane', elem_num_nodes);
    mid = N*el_coords;
    if(hit(1:7) == 'hit_off')
        line(c(:, 1),c(:, 2),c(:, 3),...
            'LineStyle','-', ...
            'LineWidth', linewidth, ...
            'Marker','none', ...
            'Color', [0.8 0.8 0.8], ...
            'MarkerSize',25, ...
            'EraseMode','none',...
            'HitTest', 'off');
    else
        line(c(:, 1),c(:, 2),c(:, 3),...
            'LineStyle','-', ...
            'LineWidth', linewidth, ...
            'Marker','none', ...
            'Color', [0.8 0.8 0.8], ...
            'MarkerSize',25, ...
            'UserData', [2 elem_num],...
            'ButtonDownFcn','perifem(''hitobject'')');
    end
    if(draw_elem_nums == 1)
        text(mid(1), mid(2), mid(3), [' ' num2str(elem_num)], 'color', 'g',...
            'UserData', [2 elem_num], 'ButtonDownFcn','perifem(''hitobject'')');
    end


   if (view_max_st == 1) 
    if(~isempty(arrows))
        center_x(1) = arrows{elem_num}(1);
        center_y(1) = arrows{elem_num}(2);
        u_max(1) = arrows{elem_num}(3);
        u_min(1) = arrows{elem_num}(4);
        v_max(1) = arrows{elem_num}(5);
        v_min(1) = arrows{elem_num}(6);

        quiver(center_x,center_y,u_max,v_max, 5000*scale_factor, 'color', 'blue')

        %quiver(center_x,center_y,u_min,v_min, 0, 'color', 'blue')
    end
   end


elseif(elem_type == 4 & (draw_elems(elem_type)==1))
    elem_num_nodes = elem_props{elem_num}(2);
    linewidth = 1;
    for node = 1 : elem_num_nodes
        el_coords(node, :) = coords(elem_con{elem_num}(node), :);
    end
    for(node = 1 : elem_num_nodes)
        if(elem_num_nodes == 4)
            if(node == 1)
                xi = [1-shrink_fact, shrink_fact/3, shrink_fact/3];
            elseif(node == 2)
                xi = [shrink_fact/3, 1-shrink_fact, shrink_fact/3];
            elseif(node == 3)
                xi = [shrink_fact/3, shrink_fact/3, 1-shrink_fact];
            elseif(node == 4)
                xi = [shrink_fact/3, shrink_fact/3, shrink_fact/3];
            end
        elseif(elem_num_nodes == 10)
            if(node == 1)
                xi = [1-shrink_fact, shrink_fact/3, shrink_fact/3];
            elseif(node == 2)
                xi = [.5, .5, shrink_fact/3];
            elseif(node == 3)
                xi = [shrink_fact/3, 1-shrink_fact, shrink_fact/3];
            elseif(node == 4)
                xi = [shrink_fact/3, .5, .5];
            elseif(node == 5)
                xi = [shrink_fact/3, shrink_fact/3, 1-shrink_fact];
            elseif(node == 6)
                xi = [.5, shrink_fact/3, .5];
            elseif(node == 7)
                xi = [.5, shrink_fact/3, shrink_fact/3];
            elseif(node == 8)
                xi = [shrink_fact/3, .5, shrink_fact/3];
            elseif(node == 9)
                xi = [shrink_fact/3, shrink_fact/3, .5];
            elseif(node == 10)
                xi = [shrink_fact/3, shrink_fact/3, shrink_fact/3];
            end
        elseif(elem_num_nodes == 8)
            if(node == 1)
                xi = [-1+shrink_fact, -1+shrink_fact, -1+shrink_fact];
            elseif(node == 2)
                xi = [1-shrink_fact, -1+shrink_fact, -1+shrink_fact];
            elseif(node == 3)
                xi = [1-shrink_fact, 1-shrink_fact, -1+shrink_fact];
            elseif(node == 4)
                xi = [-1+shrink_fact, 1-shrink_fact, -1+shrink_fact];
            elseif(node == 5)
                xi = [-1+shrink_fact, -1+shrink_fact, 1-shrink_fact];
            elseif(node == 6)
                xi = [1-shrink_fact, -1+shrink_fact, 1-shrink_fact];
            elseif(node == 7)
                xi = [1-shrink_fact, 1-shrink_fact, 1-shrink_fact];
            elseif(node == 8)
                xi = [-1+shrink_fact, 1-shrink_fact, 1-shrink_fact];
            end
        elseif(elem_num_nodes == 20)
            if(node == 1)
                xi = [-1+shrink_fact, -1+shrink_fact, -1+shrink_fact];
            elseif(node == 2)
                xi = [0, -1+shrink_fact, -1+shrink_fact];
            elseif(node == 3)
                xi = [1-shrink_fact, -1+shrink_fact, -1+shrink_fact];
            elseif(node == 4)
                xi = [1-shrink_fact, 0, -1+shrink_fact];
            elseif(node == 5)
                xi = [1-shrink_fact, 1-shrink_fact, -1+shrink_fact];
            elseif(node == 6)
                xi = [0, 1-shrink_fact, -1+shrink_fact];
            elseif(node == 7)
                xi = [-1+shrink_fact, 1-shrink_fact, -1+shrink_fact];
            elseif(node == 8)
                xi = [-1+shrink_fact, 0, -1+shrink_fact];
            elseif(node == 9)
                xi = [-1+shrink_fact, -1+shrink_fact, 0];
            elseif(node == 10)
                xi = [1-shrink_fact, -1+shrink_fact, 0];
            elseif(node == 11)
                xi = [1-shrink_fact, 1-shrink_fact, 0];
            elseif(node == 12)
                xi = [-1+shrink_fact, 1-shrink_fact, 0];
            elseif(node == 13)
                xi = [-1+shrink_fact, -1+shrink_fact, 1];
            elseif(node == 14)
                xi = [0, -1+shrink_fact, 1];
            elseif(node == 15)
                xi = [1-shrink_fact, -1+shrink_fact, 1-shrink_fact];
            elseif(node == 16)
                xi = [1-shrink_fact, 0, 1-shrink_fact];
            elseif(node == 17)
                xi = [1-shrink_fact, 1-shrink_fact, 1-shrink_fact];
            elseif(node == 18)
                xi = [0, 1-shrink_fact, 1-shrink_fact];
            elseif(node == 19)
                xi = [-1+shrink_fact, 1-shrink_fact, -1+shrink_fact];
            elseif(node == 20)
                xi = [-1+shrink_fact, 0, -1+shrink_fact];
            end
        end
        N = shape_functions(xi, 'solid', elem_num_nodes);
        c(node, :) = N*el_coords;
    end
    if(elem_num_nodes == 4)
        edges = [c(1,:); c(2,:); c(3,:); c(1,:); c(4,:); c(2,:); c(3,:); c(4,:)];
        xi = [0.25 0.25 0.25];
    elseif(elem_num_nodes == 10)
        xi = [0.25 0.25 0.25];
    elseif(elem_num_nodes == 8)
        edges = [c(1,:); c(2,:); c(3,:); c(4,:); c(1,:);...
            c(5,:); c(6,:); c(7,:); c(8,:); c(5,:);...
            c(6,:); c(2,:); c(3,:); c(7,:); c(8,:); c(4,:)];
        xi = [0 0 0];
    elseif(elem_num_nodes == 20)
        xi = [0 0 0];
    end
    N = shape_functions(xi, 'solid', elem_num_nodes);
    mid = N*el_coords;
    if(hit(1:7) == 'hit_off')
        line(edges(:, 1),edges(:, 2),edges(:, 3),...
            'LineStyle','-', ...
            'LineWidth', linewidth, ...
            'Marker','none', ...
            'Color', color, ...
            'MarkerSize',25, ...
            'EraseMode','none',...
            'HitTest', 'off');
    else
        line(edges(:, 1),edges(:, 2),edges(:, 3),...
            'LineStyle','-', ...
            'LineWidth', linewidth, ...
            'Marker','none', ...
            'Color', color, ...
            'MarkerSize',25, ...
            'EraseMode','none',...
            'UserData', [2 elem_num],...
            'ButtonDownFcn','perfifem(''hitobject'')');
    end
    if(draw_elem_nums == 1)
        text(mid(1), mid(2), mid(3), [' ' num2str(elem_num)], 'color', 'g',...
            'UserData', [2 elem_num], 'ButtonDownFcn','perifem(''hitobject'')');
    end

end

return