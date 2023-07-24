function [H, n, m, z] = METqc2sparse(qcHFileName)
    fid = fopen(qcHFileName, "r");
    n = fscanf(fid, "%d", [1 1]);
    m = fscanf(fid, "%d", [1 1]);
    z = fscanf(fid, "%d", [1 1]);
    I = sparse(eye(z));
    Z = sparse(zeros(z));
    H = sparse([]);

    for i = 1:m
        lH = sparse([]);
        for j = 1:n
            shiftList = [];
            shift = fscanf(fid, "%d", [1 1]);
            
            % Recursively read shift2 values until a empty is encountered
            while true
                shift2 = fscanf(fid, "&%d", [1 1]);
                if  isempty(shift2)
                    break;
                end
                shiftList = [shiftList, shift2];
            end
            
            if shift == -1
                lH = [lH Z];
            else
                lH = [lH, getCircshiftSum(I, shift, shiftList)];
            end
        end
        H = [H; lH];
    end

    fclose(fid);
end

function lH = getCircshiftSum(I, shift, shiftList)
    lH = circshift(I, [0 shift]);
    for k = 1:length(shiftList)
        lH = lH + circshift(I, [0 shiftList(k)]);
    end
end