function [v, wn] = myklestad_free_free(mi, EIi, dx)
    % Myklestad Method
    %--------------------------------------------------------------------------
    % Iterates through lumped masses over stations and their nearby fields to
    % find out station vectors, considering a structure free  in both edges as
    % boundary condition.
    % Returns 1) a cell array with three vectors (for the first three natural
    % frequencies), each vector containing a list of station vectors (one for
    % every station); and 2) a vector containing positive, sorted from lowest to
    % highest natural frequencies of the system.
    %
    % Input
    % ----------
    %       [mi] :      Lumped Masses Vector                              [n,1]
    %       [EIi] :     Lumped Stifness Vector                          [n+1,1]
    %       [dx]:       Field Length Vector                             [n+1,1]
    %
    % Output
    % ----------
    %       [v]:        Station Vectors                        cell array [n,3]
    %                    - three vectors (for the first three natural
    %                    frequencies), each containing station vectors for that
    %                    given frequency
    %                    - each station vector has the form [Y, psi, M, Q],
    %                    where "Y" stands for translational displacement, "psi"
    %                    for angular displacement, "M" bending moment and "Q"
    %                    shearing force
    %
    %       [wn]:       Natural Frequencies Vector                        [n,1]
    %
    % Notes
    % ----------
    % For a free-free case, we have something like
    % STATION1--FIELD1--STATION2--FIELD2--STATION3
    % In this case, we have n stations and n-1 fields
    %
    % 2) The code developed here followed the steps described in the book
    % "Fundamentals of Vibrations", by Leonard Meirovitch. The book is avaliable
    % at: http://www.iust.ac.ir/files/fnst/ssadeghzadeh_52bb7/files/EB__Fundamental_of_Vibration.pdf

    % Declaring variables used throughout this function
    num_of_stations = length(mi);
    syms w; % Symbolic variable representing a natural frequency

    % TODO: prealocate the cells with 4x4 matrices to improve perfomance
    % We need to declare this variables as cell arrays because
    % this is the only way to store matrices in a MATLAB data structure
    % See: https://www.mathworks.com/matlabcentral/answers/496101-matrix-of-matrices-matrix-with-matrices-inside-it
    TF = {}; % Field Transfer Matrices
    TS = {}; % Station Transfer Matrices
    T = {}; % Transfer Matrices

    % Remember, for a free-free case, we have something like:
    % STATION1--FIELD1--STATION2--FIELD2--STATION3
    % For STATION1, we calculate the transfer matrix as T{1} = TF{1}*T{1}
    % For STATION2, we calculate the transfer matrix as T{2} = TF{2}*T{2}
    % We don't calculate the transfer matrix at the last station (T{3}), because
    % we don't have TF{3} - to later we calculate the overall matrix and the last
    % station vector with just TS{3}
    for i = 1:(num_of_stations - 1)
        % Flexibility influence coefficients
        a_YM = (dx(i)^2) / (2 * EIi(i)); % displacement at i+1 due to a unit moment applied at i+1
        a_YQ = (dx(i)^3) / (3 * EIi(i)); % displacement at i+1 due to a unit force applied at i+1
        a_psiM = dx(i) / EIi(i); % slope at i+1 due to a unit moment applied at i+1
        a_psiQ = (dx(i)^2) / (2 * EIi(i)); % slope at i+1 due to a unit force applied at i+1

        % Transfer Matrix at nearby field
        TF{i} = [1 dx(i) a_YM -a_YQ / 2;
            0 1 a_psiM -a_psiQ;
            0 0 1 -dx(i);
            0 0 0 1];

        % Transfer Matrix at station
        TS{i} = [1 0 0 0;
            0 1 0 0;
            0 0 1 0;
            -w^2 * mi(i) 0 0 1];

        % Transfer matrix relating station vector on left side of station i+1
        % to station vector on the left side of station i
        T{i} = TF{i} * TS{i};
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATING OVERALL TRANSFER MATRIX (T)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Overall Transfer Matrix for free-free case:
    % (T_overall = TS{n}*T{n-1}*...T{2}*T{1})
    % Unfortunately, in math books the overall transfer matrix is also
    % represented by the letter "T". Here, we are a bit more explicit naming it
    % "T_overall"
    T_overall = T{1};

    for i = 2:(num_of_stations - 1)
        T_overall = T{i} * T_overall;
    end

    % TODO: abstract method to calculate a station matrix
    % so we don't repeat this here and inside the loop above
    % Final Station Transfer Matrix (TS{n})
    TS{num_of_stations} = [1 0 0 0;
                        0 1 0 0;
                        0 0 1 0;
                        -w^2 * mi(num_of_stations) 0 0 1];

    T_overall = TS{num_of_stations} * T_overall;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATING NATURAL FREQUENCIES (wn)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % This is found by carefull analysis of vRn = T*vL0
    % (where vRn is the station vector at the end edge
    % and vL0 is the station vector at the start edge)
    % Knowing that, for a free edge, Y≠0, psi≠0, M=0, and Q=0
    % we find that the determinant below must be equal to 0
    symbolic_wn = solve(det([T_overall(3, 1) T_overall(3, 2);
                        T_overall(4, 1) T_overall(4, 2)]) == 0, w);

    % Converting the natural frequencies from symbolic to numeric values
    all_wn = double(subs(symbolic_wn)); % This can contain negative numbers
    wn = sort(all_wn(all_wn > 0)); % Keep only positive numbers, sorted from lowest to highest

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % CALCULATING STATE VECTORS (v) FOR MULTIPLE MODES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % TODO: this can be a parameter of the function
    % or even selecting exactly what modes the user wants
    max_num_of_modes = 3;
    v = {}; % TODO: prealocate the cells with arrays to improve performance

    for i = 1:max_num_of_modes
        % Replacing the symbolic natural frequency with the calculated numeric
        % natural frequency for the overall matrix
        T_overall_for_wn = double(subs(T_overall, w, wn(i)));

        % We'll have one station vector for every station:
        % If we have: STATION1--FIELD1--STATION2--FIELD2--STATION3
        % We get:     v1----------------v2----------------v3

        % Preallocating the variable that will contain the station vectors
        % with 4x1 vectors
        v{i} = {};

        for j = 1:(num_of_stations)
            v{i}{j} = zeros(4, 1);
        end

        % For a free edge:
        M0 = 0; Q0 = 0;

        % We can arbitrate any value as the initial value of psi0, because the scale
        % of the station vector doesn't matter
        psi0 = 1;

        % The formula to determine Y0 is also found by analysis of the boundary
        % conditions in the equation vRn = T_overall*vL0.
        Y0 = -T_overall_for_wn(3, 2) * psi0 / T_overall_for_wn(3, 1);

        v{i}{1} = [Y0; psi0; M0; Q0]; % v1, station vector at start/left edge (wall)

        % Then:
        % v{i}{2} = T{1}*v{1}, equivalent to v2 = T1*v1
        % v{i}{3} = T{2}*v{2}, equivalent to v3 = T2*v2
        % ...
        % v{i}{n} = T{n-1}*v{n-1}, equivalent to vn = Tn*vn
        % which is also equivalent to vn = T*v0
        % (vRn = T_overall*vL0 mentioned before)
        for j = 2:(length(v{i}))
            % Replacing the symbolic natural frequency with the calculated numeric
            % natural frequency for the transfer matrix at given station
            T_numeric = double(subs(T{j - 1}, w, wn(i)));
            v{i}{j} = T_numeric * v{i}{j - 1};
        end

    end

end
