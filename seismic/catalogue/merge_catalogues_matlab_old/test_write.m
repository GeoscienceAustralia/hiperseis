for i=22015:length(ev(:,1))
    % Write the event info
    %                        1       2       3        4      5       6        7     8       9       10     11     12     13     14     15     16
    %                       %YY      MM      DD      hh      mm      ss      lo     la      de      ph     mb,    ms     ml,    mw,    ID     GAP    EVnumber
    if ev(i,17) == 1
        fprintf(fid,'#EHB, %04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %06.3f, %10.5f, %10.5f, %10.5f, %5.0f, %4.2f, %4.2f, %4.2f, %4.2f, %7.0f, %7.4f, %1.0f\n',...
            ev(i,1:16), i );
    elseif ev(i,17) == 2
        fprintf(fid,'#ISC, %04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %06.3f, %10.5f, %10.5f, %10.5f, %5.0f, %4.2f, %4.2f, %4.2f, %4.2f, %7.0f, %7.4f, %1.0f\n',...
            ev(i,1:16), i );
    elseif ev(i,17) == 3
        fprintf(fid,'#USGS, %04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %06.3f, %10.5f, %10.5f, %10.5f, %5.0f, %4.2f, %4.2f, %4.2f, %4.2f, %7.0f, %7.4f, %1.0f\n',...
            ev(i,1:16), i );
    elseif ev(i,17) == 4
        fprintf(fid,'#GG, %04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %06.3f, %10.5f, %10.5f, %10.5f, %5.0f, %4.2f, %4.2f, %4.2f, %4.2f, %7.0f, %7.4f, %1.0f\n',...
            ev(i,1:16), i );
    elseif ev(i,17) == 5
        fprintf(fid,'#GA, %04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %06.3f, %10.5f, %10.5f, %10.5f, %5.0f, %4.2f, %4.2f, %4.2f, %4.2f, %7.0f, %7.4f, %1.0f\n',...
            ev(i,1:16), i );
    end
    l = l+1; % line number
    mod1000 = mod(l,1000);
    if mod1000 == 1
        disp(['Event number: ',num2str(i)])
    end
    if isempty(picks(i).at)
        
    else
        for j=1:length(picks(i).at(:,1))
            fprintf(fid,'%s, ',picks(i).st(j,:)); % 1. stations
            if isempty(picks(i).ch)
                fprintf(fid,'%s, ',' ');
            else
                fprintf(fid,'%s, ',picks(i).ch(j,:)); % 2. channel
            end
            if isempty(picks(i).lc)
                fprintf(fid,'%s, ',' ');
            else
                fprintf(fid,'%s, ',picks(i).lc(j,:)); % 3. location code 00, 10, ...
            end
            if isempty(picks(i).ne)
                fprintf(fid,'%s, ',' ');
            else
                fprintf(fid,'%s, ',picks(i).ne(j,:)); % 4. network
            end
            if isempty(picks(i).lo)
                fprintf(fid,'%s, ',' ');
            else
                fprintf(fid,'%f, ',picks(i).lo(j)); % 5. lon
            end
            if isempty(picks(i).la)
                fprintf(fid,'%s, ',' ');
            else
                fprintf(fid,'%f, ',picks(i).la(j)); % 6. lat
            end
            if isempty(picks(i).el)
                fprintf(fid,'%s, ',' ');
            else
                fprintf(fid,'%f, ',picks(i).el(j)); % 7. ele
            end
            fprintf(fid,'%s, ',picks(i).ph(j,:)); % 8. phase
            fprintf(fid,'%04.0f, %02.0f, %02.0f, %02.0f, %02.0f, %05.2f, ',picks(i).at(j,:)); % 9. 10. 11. 12. 13. 14. Arrival Time
            fprintf(fid,'%7.3f \n',picks(i).dist(j)); % 15. distance
            l=l+1;
        end
    end
end