function Output = SimmonsModel (OutType, Input, WF_Au_eV)
    % Simmons' Model (tunneling through a trapezoidal barrier)

    %Constants:
    %http://physics.nist.gov/cuu/Constants/index.html
    % h-bar = 1.054 571 726(47) x 10-34 Js  [2013-12] Planck constant over 2 pi (= reduced Planck constant)
    % m_e = 9.109 382 91(40) x 10-31 kg     [2013-12] Electron mass
    % eV = 1.602 176 565(35) x 10-19 J      [2013-12] Electron volt

    %Work function for Au (in J, not eV!)
    % http://en.wikipedia.org/wiki/Work_function
    % "The work function of gold", Sachtler et al.: http://dx.doi.org/10.1016/0039-6028(66)90083-5
    WF_Au = WF_Au_eV * 1.60217656535e-19;

    %Semilog plot:
    % figure(),semilogy( linspace(0,5e-10,100), Simmons_Transmission(linspace(0,5e-10,100)) )

    switch OutType
        case 'Transmission'
            %Calculate transmission probability:
            x = Input;
            Output = exp( -2/1.05457172647e-34 * sqrt(2*9.1093829140e-31) * sqrt(WF_Au) * x );
        case 'Distance'
            %Calculate potential barrier width (= electrode distance):
            T = Input;
            Output = log(T)/( -2/1.05457172647e-34 * sqrt(2*9.1093829140e-31) * sqrt(WF_Au) );
        otherwise
            error('Unknown OutType for SimmonsModel function!');
    end
end
