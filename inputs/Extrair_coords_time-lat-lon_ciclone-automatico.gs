* ============================================================================
* Autor do código:
*  ____   ___  _   _    _    _     ____  
* |  _ \ / _ \| \ | |  / \  | |   |  _ \ 
* | |_) | | | |  \| | / _ \ | |   | | | |
* |  _ <| |_| | |\  |/ ___ \| |___| |_| |
* |_| \_\\___/|_| \_/_/   \_\_____|____/ 
*
*        >>>  R O N A L D   G U I U S E P P I   R A M Í R E Z   N I N A  <<<
*              Ph.D. Candidate in Atmospheric Sciences
*
* Disciplina: Meteorologia Sinótica III
*
* Professores:
*    - Prof. Dr. Pedro Leite da Silva Dias
*    - Prof. Dr. Ricardo Hallak
*
* Universidade de São Paulo (USP)
* Instituto de Astronomia, Geofísica e Ciências Atmosféricas (IAG-USP)
* Departamento de Ciências Atmosféricas
*
* ============================================================================
* EXTRAI CENTRO DO CICLONE
* ============================================================================
* Reiniciar o terminal se existir algum processo passado 
* Configurar a tela de cor branca
'reinit'
'set display color white'
'c'

* ============================================================================
* Abre o arquivo
* ============================================================================
'sdfopen mslp_ciclone_dec_2025.nc'

* ============================================================================
* Loop no tempo (n) - ALTERAR O N DE ACORDO COM SEU NUMERO DE TEMPOS DE ANALISES
* ============================================================================
nt = 72
n = 1

* ============================================================================
* Área de busca do ciclone
* Ajuste conforme necessário
* ============================================================================
lat1 = -45
lat2 = -10
lon1 = -65
lon2 = -25
* ============================================================================
* Arquivo de saída
* ============================================================================
outfile = 'track_ciclone_2025_automatico.txt'
*rc = write(outfile, 'data latitude longitude pmin_hPa')
rc = write(outfile, 'time;Lat;Lon')

* ============================================================================
* Loop temporal
* ============================================================================
while (n <= nt)

'set t ' n
'c'

* ============================================================================
* Definir o domínio de visualização
* ============================================================================
'set lat ' lat1 ' ' lat2
'set lon ' lon1 ' ' lon2
*'set mpdset hires'
'set mpdset hires brmap'

* Obter a data
'q time'
 data = subwrd(result,3)
    
* ==========================================
* Transformar o formato do tempo
* ==========================================
tstr = subwrd(result,3)
* ============================================================================
* Extrair componentes
* Exemplo: 09Z05DEC2025
* ============================================================================
hh  = substr(tstr,1,2)
dd  = substr(tstr,4,2)
mon = substr(tstr,6,3)
yyyy = substr(tstr,9,4)
* ============================================================================
* Converter mês para número
* ============================================================================
if (mon = 'JAN'); mm = '01'; endif
if (mon = 'FEB'); mm = '02'; endif
if (mon = 'MAR'); mm = '03'; endif
if (mon = 'APR'); mm = '04'; endif
if (mon = 'MAY'); mm = '05'; endif
if (mon = 'JUN'); mm = '06'; endif
if (mon = 'JUL'); mm = '07'; endif
if (mon = 'AUG'); mm = '08'; endif
if (mon = 'SEP'); mm = '09'; endif
if (mon = 'OCT'); mm = '10'; endif
if (mon = 'NOV'); mm = '11'; endif
if (mon = 'DEC'); mm = '12'; endif
* ============================================================================
* Montar formato final
* YYYY-MM-DD-HHmm
* ============================================================================
data_fmt = yyyy '-' mm '-' dd '-' hh '00'
* ============================================ 
    
* Fazer um plot do campo de pressão reduzida ao nível do mar 
'set gxout contour'
'set cint 2'
'd msl/100'
'draw title PRNMM ' data

* ============================================================================
* Valor mínimo da pressão na área
* ============================================================================
'd amin(msl,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
pmin = subwrd(result,4)/100.

* ============================================================================
* Longitude do mínimo
* ============================================================================
'd aminlocx(msl,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
xlon = subwrd(result,4)
'set x ' xlon
longitude = subwrd(result,4)

* ============================================================================
* Latitude do mínimo
* ============================================================================
'd aminlocy(msl,lon='lon1',lon='lon2',lat='lat1',lat='lat2')'
ylat = subwrd(result,4)
'set y ' ylat
latitude = subwrd(result,4)

* ============================================================================
* Imprimir na tela as coordenadas geográficas
* ============================================================================
say '----------------------------------------'
say 'Tempo: ' data
say 'Centro minimo: lat=' latitude ' lon=' longitude ' p=' pmin ' hPa'

* ============================================================================
* Salvar o arquivo .txt com a informação para o track
* ============================================================================
*linha = data ' ' latitude ' ' longitude ' ' pmin
linha = data_fmt ';' latitude ';' longitude
rc = write(outfile, linha)

n = n + 1

endwhile

rc = close(outfile)

say ' '
say 'Processamento finalizado.'
say 'Saida salva em: ' outfile
