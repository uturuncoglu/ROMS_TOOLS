;***********************************************************
;*** Ufuk Turuncoglu ***************************************
;*** Send bug reports to turuncu@be.itu.edu.tr *************
;***********************************************************
;-----------------------------------------------------------
; List of available functions/procedures 
; in alphabetical order
;-----------------------------------------------------------
; ROMS_u2rho - converts u point data to rho
; ROMS_v2rho - converts v point data to rho
;-----------------------------------------------------------

undef("ROMS_u2rho")
function ROMS_u2rho(input:numeric)
local output, dims, rank
begin
  ;--- get variable dimension and rank ---
  dims = dimsizes(input)
  rank = dimsizes(dims)

  ;--- move data from rho to u points ---
  if (rank .eq. 4) then
    Mp  = dims(2)
    L = dims(3)
    Lp = L+1
    Lm = L-1
    output = new((/ dims(0), dims(1), Mp, Lp /), typeof(input), input@_FillValue)

    output(:,:,:,1:L-1) = 0.5*(input(:,:,:,0:Lm-1)+input(:,:,:,1:L-1))
    output(:,:,:,0) = output(:,:,:,1)
    output(:,:,:,Lp-1) = output(:,:,:,L-1)
  end if

  if (rank .eq. 3) then
    Mp  = dims(1)
    L = dims(2)
    Lp = L+1
    Lm = L-1
    output = new((/ dims(0), Mp, Lp /), typeof(input), input@_FillValue)

    output(:,:,1:L-1) = 0.5*(input(:,:,0:Lm-1)+input(:,:,1:L-1))
    output(:,:,0) = output(:,:,1)
    output(:,:,Lp-1) = output(:,:,L-1)
  end if

  if (rank .eq. 2) then
    Mp  = dims(0)
    L = dims(1)
    Lp = L+1
    Lm = L-1
    output = new((/ Mp, Lp /), typeof(input), input@_FillValue)

    output(:,1:L-1) = 0.5*(input(:,0:Lm-1)+input(:,1:L-1))
    output(:,0) = output(:,1)
    output(:,Lp-1) = output(:,L-1)
  end if

  ;--- copy coordinate attribute --
  copy_VarMeta(input, output)

  return(output)
end

undef("ROMS_v2rho")
function ROMS_v2rho(input:numeric)
local output, dims, rank
begin
  ;--- get variable dimension and rank ---
  dims = dimsizes(input)
  rank = dimsizes(dims)

  ;--- move data from rho to u points ---
  if (rank .eq. 4) then
    M  = dims(2)
    Mp = M+1
    Mm = M-1
    Lp = dims(3)
    output = new((/ dims(0), dims(1), Mp, Lp /), typeof(input), input@_FillValue)

    output(:,:,1:M-1,:) = 0.5*(input(:,:,0:Mm-1,:)+input(:,:,1:M-1,:))
    output(:,:,0,:) = output(:,:,1,:)
    output(:,:,Mp-1,:) = output(:,:,M-1,:)
  end if

  if (rank .eq. 3) then
    M  = dims(1)
    Mp = M+1
    Mm = M-1
    Lp = dims(2)
    output = new((/ dims(0), Mp, Lp /), typeof(input), input@_FillValue)

    output(:,1:M-1,:) = 0.5*(input(:,0:Mm-1,:)+input(:,1:M-1,:))
    output(:,0,:) = output(:,1,:)
    output(:,Mp-1,:) = output(:,M-1,:)
  end if

  if (rank .eq. 2) then
    M  = dims(0)
    Mp = M+1
    Mm = M-1
    Lp = dims(1)
    output = new((/ Mp, Lp /), typeof(input), input@_FillValue)

    output(1:M-1,:) = 0.5*(input(0:Mm-1,:)+input(1:M-1,:))
    output(0,:) = output(1,:)
    output(Mp-1,:) = output(M-1,:)
  end if

  ;--- copy coordinate attribute --
  copy_VarMeta(input, output)

  return(output)
end
