function switch_case(cluster,n,B)

    switch cluster
        case {3}
            for i=1:1:n
                if B(i)~=0
                    plot(x(i),y(i),'yp');
                    hold on;
                end
            end
        
    end
    
    