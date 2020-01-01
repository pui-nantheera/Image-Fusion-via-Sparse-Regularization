function q=q_abf(a,b,f,p)

% Objective Fusion Performance measure for fusion of 2 multisensor images (A and B)
% into a fused image F (all images must be of the same size)
% returns a single value for q
% p is a set of 4 parameters for the sigmoid non-linearities,
% inputing p=[] will use default parameters
%     q=q_abf(a,b,f,p)
% VP, October 2000, revision Jul 2003.

% Ucitavamo slike Reading images
if ischar(a), a=double(imread(a)); end
if ischar(b), b=double(imread(b)); end
if ischar(f), f=double(imread(f)); end

if isempty(p)
  p=[0.8 24 0.7 11];     % Standardni parametri: Standard parameters
end

% Odredjujemo parametre nelinearnosti: Calculating non-linear parameters
sa=p(1);
ka=p(2);
sg=p(3);
kg=p(4);
gg=(1+exp(-kg*(1-sg)));
ga=(1+exp(-ka*(1-sa)));

if length(size(a))>2
    a=(a(:,:,1)+a(:,:,2)+a(:,:,3))/3;
end
if length(size(b))>2
    b=(b(:,:,1)+b(:,:,2)+b(:,:,3))/3;
end
if length(size(f))>2
    f=(f(:,:,1)+f(:,:,2)+f(:,:,3))/3;
end

% Parametri slika : Image Parameters
[visina sirina]=size(a);
% Racunamo parametre ivicnih elemnenata za sve slike: 
% We are computing parameters of edge elements for all images
[a_j,a_o]=i_parametri(a);
[b_j,b_o]=i_parametri(b);
[f_j,f_o]=i_parametri(f);



% Odredjujemo prag detekcije i korigujemo vrednosti jacine: 
% Calculating detection Threshold and correcting the edge strength
dt=2;
gmin=1;

prag_a=a_j<dt;          % Prag detekcije: Detection Threshold
prag_b=b_j<dt;
prag_f=f_j<dt;
a_j=a_j.*(1-prag_a)+gmin*prag_a;  %***** a_j(prag_a) = gmin; as simple as that
b_j=b_j.*(1-prag_b)+gmin*prag_b;
f_j=f_j.*(1-prag_f)+gmin*prag_f;

a_o=a_o-(prag_a.*prag_f).*(a_o-f_o); %***** this can be removed if done as coment in line 71
b_o=b_o-(prag_b.*prag_f).*(b_o-f_o); 

% Ulaz u sistem mere (OEFP): Entering the measurement system (OEFP):
% Racunamo linearne vrednosti opadanja/rasta jacine i orijentacije :
% Calculating the linear values of desrease/increase of the strength and
% orientation of the edge values in A and B in F
% ivicnih vrednosti iz A i B u F
marg=0.000034;
g_af=(f_j+marg)./(a_j+marg);
g_bf=(f_j+marg)./(b_j+marg);
a_af=abs(abs(a_o-f_o)-pi/2)/(pi/2);  % after this add:  a_af(prag_a.*prag_f) = 1;
a_bf=abs(abs(b_o-f_o)-pi/2)/(pi/2);

prag_gaf=g_af>1;
prag_gbf=g_bf>1;
prag_sv=prag_gaf.*prag_gbf.*prag_a.*prag_b;

g_af=g_af.*(1-prag_gaf)+prag_gaf./g_af;  %***** g_af = (min(f_j,a_j)+eps)/(max(f_j,a_j)+eps);
g_bf=g_bf.*(1-prag_gbf)+prag_gbf./g_bf;
a_j=a_j.*(1-prag_sv)+f_j.*prag_sv;  %**** why? specially why after having use a_j already...
b_j=b_j.*(1-prag_sv)+f_j.*prag_sv;  %**** why?

clear prag_a prag_b prag_f prag_gaf prag_gbf prag_sv

% Propustamo faktore odrzanja jacine i orijentacije kroz nelinearnosti
% We are filtering the strength and orientation through non-linearities
qg_af=gg./(1+exp(-kg*(g_af-sg)));
qg_bf=gg./(1+exp(-kg*(g_bf-sg)));

qa_af=ga./(1+exp(-ka*(a_af-sa)));
qa_bf=ga./(1+exp(-ka*(a_bf-sa)));

% Racunamo koeficijent odrzanja informacija za celokupne
% ivicne elemente Qaf i Qbf: Calculating coeficient of retained information
% for all edge elements  in Qaf and Qbf
q_af=sqrt(qg_af.*qa_af);
q_bf=sqrt(qg_bf.*qa_bf);

% Racunamo provizorne mape znacaja: Calculating provisional maps of
% importance
wa=a_j;
wb=b_j;              % Racunamo konacne mape znacaja: Calc. the final maps of importance
W=sum(sum(wa+wb));

clear za zb qg_af qb_bf qa_af qa_bf a_j b_j f_j a_o b_o f_o

% Racunamo finalnu meru sjedinjavanja: Calculating the final measure of the
% fusion
q=sum(sum(q_af.*wa+q_bf.*wb))/W;
