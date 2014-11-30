;-----------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "ROMS_func.ncl"
;-----------------------------------------------------------
;***********************************************************
;*** Ufuk Turuncoglu ***************************************
;*** Send bug reports to turuncu@be.itu.edu.tr *************
;***********************************************************
begin
  ;--- parameters ---
  ifile = "ocean_avg.nc"
  ofile = "ocean_avg_zlev.nc"
  vlist = (/ "u"      , "v"         , \
            "lon_rho" , "lat_rho"   , \
            "mask_rho", "ocean_time" /)
  nvar = dimsizes(vlist)
  zlev = (/ 0, 2.5, 5.0, 10.0, 15.0 /)
  nlev = dimsizes(zlev)

  ;---------------------------------------------------------
  ; ********************************************************
  ; *** Main section ***************************************
  ; ********************************************************
  ;
  ; Do not change the code after this point if it is not 
  ; really necessary !!! 
  ;---------------------------------------------------------

  ;--- open input dataset ---
  fin = addfile(ifile, "r")

  ;--- read grid related variables ---
  s_rho = fin->s_rho
  s_w = fin->s_w
  cs_r = fin->Cs_r
  cs_w = fin->Cs_w
  vt = fin->Vtransform
  vs = fin->Vstretching
  h = fin->h
  hc = fin->hc

  msk2d = fin->mask_rho
  lat2d = fin->lat_rho
  lon2d = fin->lon_rho   
  time = fin->ocean_time

  ;--- get dimensions ---
  dims = dimsizes(msk2d)
  nlat = dims(0)
  nlon = dims(1)    
  km = dimsizes(s_rho)
  ntime = dimsizes(time)

  ;---------------------------------------------------------
  ; STEP 1: Create output file
  ;---------------------------------------------------------

  system("/bin/rm -f "+ofile)
  ncout = addfile(ofile, "c")

  ;--- define mode on ---
  setfileoption(ncout, "DefineMode", True)

  ;--- define dimensions ---
  do i = 0, nvar-1
    dimNames = getfilevardims(fin, vlist(i))
    dimSizes = getfilevardimsizes(fin, vlist(i))

    indx = ind(dimNames .eq. "s_rho")
    if (.not. ismissing(indx)) then
      dimNames(indx) = "lev"
      dimSizes(indx) = nlev
    end if
    delete(indx)

    indx = ind(dimNames .eq. "eta_u" .or. dimNames .eq. "eta_v")
    if (.not. ismissing(indx)) then
      dimNames(indx) = "eta_rho"
      dimSizes(indx) = nlat 
    end if
    delete(indx)

    indx = ind(dimNames .eq. "xi_u" .or. dimNames .eq. "xi_v")
    if (.not. ismissing(indx)) then
      dimNames(indx) = "xi_rho"
      dimSizes(indx) = nlon 
    end if
    delete(indx)

    dimUnlim = new((/ dimsizes(dimNames) /), "logical")
    dimUnlim = where(dimNames .eq. "ocean_time", True, False)
    filedimdef(ncout, dimNames, dimSizes, dimUnlim)
    delete([/ dimNames, dimSizes, dimUnlim /])
  end do

  ;--- define variables ---
  do i = 0, nvar-1
    dimNames = getfilevardims(fin, vlist(i))

    indx = ind(dimNames .eq. "s_rho")
    if (.not. ismissing(indx)) then
      dimNames(indx) = "lev"
    end if
    delete(indx)
  
    indx = ind(dimNames .eq. "eta_u" .or. dimNames .eq. "eta_v")
    if (.not. ismissing(indx)) then
      dimNames(indx) = "eta_rho"
    end if
    delete(indx)    

    indx = ind(dimNames .eq. "xi_u" .or. dimNames .eq. "xi_v")
    if (.not. ismissing(indx)) then
      dimNames(indx) = "xi_rho"
    end if
    delete(indx)
    
    varType = getfilevartypes(fin, vlist(i))
    filevardef(ncout, vlist(i), varType, dimNames)
    delete(dimNames)
  end do

  ;--- add attributes ---
  do i = 0, nvar-1
    varAtts = getfilevaratts(fin, vlist(i))
    dimNames = getfilevardims(fin, vlist(i))
    varType = getfilevartypes(fin, vlist(i))
    ndims = dimsizes(dimNames)

    if (ndims .eq. 1) then
      dumm = fin->$vlist(i)$(0)
    end if
    if (ndims .eq. 2) then
      dumm = fin->$vlist(i)$(0,0)
    end if
    if (ndims .eq. 3) then
      dumm = fin->$vlist(i)$(0,0,0)
    end if    
    if (ndims .eq. 4) then
      dumm = fin->$vlist(i)$(0,0,0,0)
    end if

    if (varType .eq. "float") then
      atts = 0.0
    end if
    if (varType .eq. "double") then
      atts = 0.0d0
    end if

    do j = 0, dimsizes(varAtts)-1
      if (varAtts(j) .eq. "coordinates") then
        str = dumm@$varAtts(j)$
        str = str_sub_str(str, "_u", "_rho") 
        str = str_sub_str(str, "_v", "_rho")
        str = str_sub_str(str, "s_rho", "lev")
        atts@$varAtts(j)$ = str
      else 
        atts@$varAtts(j)$ = dumm@$varAtts(j)$
      end if
    end do
    filevarattdef(ncout, vlist(i), atts)
    delete([/ varAtts, dimNames, varType, dumm, atts /])
  end do

  ;--- add variable for vertical levels ---
  filevardef(ncout, "lev" , "double", "lev")
  attr = True
  attr@long_name = "Depth"
  attr@units = "meters" 
  filevarattdef(ncout, "lev", attr)
  delete(attr)

  ;--- exit file definition mode ---
  setfileoption(ncout, "DefineMode", False)

  ;--- fill variables without vertical dimension---
  do i = 0, nvar-1
    dimNames = getfilevardims(fin, vlist(i))
    if (all(dimNames .ne. "s_rho")) then
      ncout->$vlist(i)$ = (/ fin->$vlist(i)$ /)
    end if
    delete(dimNames)
  end do
  ncout->lev = (/ zlev /) 

  ;---------------------------------------------------------
  ; STEP 2: Perform interpolation from s layers to z levels
  ;---------------------------------------------------------

  zr = new((/ km, nlat, nlon /), "double", 1.0d20)

  do i = 0, ntime-1
    ;--- time measurement ---
    time1 = get_cpu_time()

    ;--- read data ---
    zeta = fin->zeta(i,:,:)

    ;--- calculate depths ---
    if (vt .eq. 1) then 
      do k = 0, km-1
        z0 = hc*s_rho(k)+(h-hc)*cs_r(k)
        zr(k,:,:) = z0+zeta*(1.0d0+z0/h)
      end do
    end if
    if (vt .eq. 2) then
      do k = 0, km-1
        z0 = (hc*s_rho(k)+h*cs_r(k))/(hc+h)
        zr(k,:,:) = zeta+(zeta+h)*z0
      end do
    end if

    ;--- perform interpolation ---
    do j = 0, nvar-1
      dimNames = getfilevardims(fin, vlist(j))

      if (any(dimNames .eq. "s_rho")) then
        varType = getfilevartypes(fin, vlist(j))
        ndims = dimsizes(dimNames)
        print(j+" "+vlist(j)+" "+ndims+" "+varType)

        if (ndims .eq. 4) then 
          flag = True
          if (vlist(j) .eq. "u") then
            flag = False
            vars = ROMS_u2rho(fin->$vlist(j)$(i,:,:,:))
          end if
          if (vlist(j) .eq. "v") then
            flag = False
            vars = ROMS_v2rho(fin->$vlist(j)$(i,:,:,:))
          end if
          if (flag) then
            vars = fin->$vlist(j)$(i,:,:,:)
          end if
        end if
        if (ndims .eq. 3) then 
          flag = True
          if (vlist(j) .eq. "u") then      
            flag = False
            vars = ROMS_u2rho(fin->$vlist(j)$)
          end if
          if (vlist(j) .eq. "v") then
            flag = False
            vars = ROMS_v2rho(fin->$vlist(j)$)
          end if
          if (flag) then
            vars = fin->$vlist(j)$
          end if
        end if

        varz = (/ linint1_n_Wrap(-zr(::-1,:,:), vars(::-1,:,:), False, zlev, 0, 0) /)
        varz(0,:,:) = where(ismissing(varz(0,:,:)), vars(km-1,:,:), varz(0,:,:))

        if (ndims .eq. 3) then
          ncout->$vlist(j)$ = (/ varz /)        
        end if
        if (ndims .eq. 4) then
          ncout->$vlist(j)$(i,:,:,:) = (/ varz /)        
        end if
      end if

      delete(dimNames)         
    end do

    ;--- time measurement ---
    time2 = get_cpu_time()

    print("[debug] -- time step '"+sprinti("%08d", i)+"' processed in "+\
          sprintf("%8.2f", (time2-time1))+" sec") 
  end do
end
