# :Title: coord.rb: Geographical coordinate transformations.
# coord.rb: Geographical coordinate transformations.
# Author:: Matt Foster (mailto:mattfoster@clara.co.uk)
# Copyright:: Copyright (c) 2005 Matt Foster
#
# A Collection of classes, modules and extensions which implement geographical
# coordinate transforms.
#

# Extensions to the Numeric class for conversion from degrees to radians and 
# radians to degrees.
class Numeric
  def rad_to_deg
    self * (180/Math::PI)
  end

  def deg_to_rad
    self * (Math::PI/180)
  end
end

# Coordinate Conversions methods, for mixing into things. 
# Contains verious methods for converting single coordinates from one type to
# another.
module CoordinateConversion
  
  # This function exists so that the classes which it is being mixed into
  # knowns what 
  def listTypes
    @types = ['dms', 'dm', 'dd']
  end

  # This doesn't so anything, but might come in handy.
  def dms_to_dms(d, m, s)
    [d, m, s]
  end

  # This doesn't so anything, but might come in handy.
  def dd_to_dd(dd)
    dd
  end

  # This doesn't so anything, but might come in handy.
  def dm_to_dm(d, m)
    [d, m]
  end

  # Convers DMS to decimal degrees (D.d).
  def dms_to_dd(d, m, s)
    d + m.to_f/60 + s.to_f/3600
  end

  # Convert D.d to DMS
  def dd_to_dms(dd)
    d = dd.floor
    f = (dd - d) * 60
    m = f.floor
    f = (f - m) * 60
    s = f

    [d, m, s.to_f]
  end

  # Convert from DMS to DM.m
  def dms_to_dm(d, m, s)
    f = s/60.0
    [d, m+f]
  end

  # Convert from DM.m to D.d
  def dm_to_dd(d, m)
    f = m/60.0
    d+f
  end

  # Convert from DD to DM.m
  def dd_to_dm(d)
    f = (d - d.floor)*60
    [d.floor, f]
  end

  # Convert from DM to DMS
  def dm_to_dms(d, m)
    f = (m - m.floor)*60
    [d, m.floor, f]
  end
end

# Transforms sets of coordinates:
# Currently supports:
# * lat + long (+ ellipsoid height)
# * 3D cartesian (x, y, z)
# * Transverse Mercator (easting, northing)
class Coordinates

  
  @@m_eqn = %q{
    m_1 	= b * f_0;
    m_2	    = (1+n+(5.0/4.0)*n**2 + (5.0/4.0)*n**3) * phi_d;
    m_2     -= (3*n + 3*n**2+(21.0/8.0)*n**3) * sin(phi_d) * cos(phi_p)
    m_2     += ((15.0/8.0)*n**2+(15.0/8.0)*n**3) * sin(2*phi_d) * cos(2*phi_p)
    m_2	    -= ((35.0/24.0)*n**3) * sin(3*phi_d) * cos(3*phi_p);
    m 		= m_1 * m_2
  }

  def initialize(debug = false)
    @debug = debug
    
    # Generate a load of trig functions: just sugar for the mass of
    # calculations in ll_to_utm.
    ['sin', 'cos', 'tan'].each {|op|
      eval "def #{op}(v); Math.#{op}(v); end"
      (1..8).each {|x|
        eval "def #{op}#{x}(v); (Math.#{op}(v))**#{x}; end"
      }
    }
  end

  # Convert from lat (+phi+ (rads)) and long (+lam+ (rads)) to easting and
  # northing
  def ll_to_en(phi, lam, proj_type = 'National Grid')
   
    proj = ProjectionInfo.new(proj_type)
    phi_0, lam_0, n_0, e_0, f_0, a, b, e_2 = proj.returnInfo
    
    # A few useful differences.
    phi_d   = phi - phi_0
    phi_p   = phi + phi_0
    lam_d   = lam - lam_0
    lam_p   = lam + lam_0

    n 		= (a - b) / (a + b)

    v 		= a * f_0 * (1-e_2 * sin2(phi))**-0.5
    rho 	= a * f_0 * (1-e_2) * (1-e_2*sin2(phi))**-1.5

    eta_2   = (v / rho) - 1

    m = ''
    eval(@@m_eqn)

    i		    = m + n_0
    ii		= (v/2.0) * sin(phi) * cos(phi)
    iii		= (v/24.0)* sin(phi) * cos3(phi)
    iii     *= (5 - tan2(phi) + 9 * eta_2)
    iiia  	= (v/720.0) * sin(phi) * cos5(phi)
    iiia    *= 61 - 58 * tan2(phi) + tan4(phi)
    iv		= v* cos(phi)
    vv		= (v/6.0) * cos3(phi) * ((v/rho)- tan2(phi))

    va		= (v/120.0) * cos5(phi) 
    vb      = 5 - 18 * tan2(phi) + tan4(phi) + 14 * eta_2
    vb      -= (58 * tan2(phi) * eta_2)
    vi      = va * vb

    n		= i + ii * lam_d**2 + iii * lam_d**4 + iiia * lam_d**6
    e		= e_0 + iv * lam_d + vv * lam_d**3 + vi * lam_d**5

    if @debug
      [:v, :rho, :eta_2, :m, :i, :ii, :iii, :iiia, :iv, :vv, :e, :n].each {|v|
        name = n.id2name
        val =  eval(n.id2name)
        printf "%s\t%10.10e\n", name, val
      }
    end

    [e, n]
  end

  # Convert easting (+ee+) and northing (+nn+) to lat and long.
  def en_to_ll(ee, nn, proj_type = 'National Grid')

    proj    = ProjectionInfo.new(proj_type)
    phi_0, lam_0, n_0, e_0, f_0, a, b, e_2 = proj.returnInfo

    n = (a-b) / (a+b)
    phi_pr  = (nn - n_0) / (a * f_0) + phi_0
    
    # A few useful differences.
    phi_d   = phi_pr - phi_0
    phi_p   = phi_pr + phi_0
    
    m = 0
    eval(@@m_eqn)
    comp = nn - n_0 - m
     
    while comp >= 0.01e-3
      phi_pr = comp / (a * f_0) + phi_pr
      phi_d = phi_pr - phi_0
      phi_p = phi_pr + phi_0
      m = ''
      eval(@@m_eqn)
      comp  = nn - n_0 - m
    end

    v 	    = a * f_0 * (1-e_2 * sin2(phi_pr))**-0.5
    rho 	= a * f_0 * (1-e_2) * (1-e_2*sin2(phi_pr))**-1.5
    eta_2   = (v / rho) - 1

    vii     = tan(phi_pr) / (2*rho*v)
    viiia   = tan(phi_pr) / (24*rho*v**3)
    viiib   = 5+3*tan2(phi_pr) + eta_2 - 9*(tan2(phi_pr))*eta_2
    viii    = viiia * viiib

    ixa     = tan(phi_pr) / (720*rho*v**5) 
    ixb     = 61 + 90*tan2(phi_pr) + 45*tan4(phi_pr)
    ix      = ixa * ixb

    x       = 1 / (cos(phi_pr)*v)

    xia     = 1 / (cos(phi_pr)*6*v**3)
    xib     = v / rho + 2*tan2(phi_pr)
    xi      = xia * xib

    xiia    = 1 / (cos(phi_pr)*120*v**5)
    xiib    = 5 + 28*tan2(phi_pr)+24*tan4(phi_pr)
    xii     = xiia * xiib

    xiiia   = 1 / (cos(phi_pr)*5040*v**7)
    xiiib   = 61 + 662*tan2(phi_pr)+1320*tan4(phi_pr)+720*tan6(phi_pr)
    xiii    = xiiia * xiiib

    # Fancy debug printing.
    if @debug
      [:v, :rho, :eta_2, :vii, :viii, :ix, :x, :xi, :xii, :xiii].each {|n|
        name = n.id2name
        val =  eval(n.id2name)
        printf "%s\t%10.10e\n", name, val
      }
    end

    phi     = phi_pr - vii*(ee-e_0)**2 + viii*(ee-e_0)**4 - ix*(ee-e_0)**6
    lam     = lam_0 + x*(ee-e_0) - xi*(ee-e_0)**3 + xii*(ee-e_0)**5 - xiii*(ee-e_0)**7
    [phi, lam]
  end

  # Convert easting and northing to national grid.
  def en_to_ng(e, n, grid_type = 'GB')

    g = GridInfo.new(grid_type)
    grid = g.returnInfo
    
    x = (e / 100000).floor
    y = (n / 100000).floor

    new_easting = (e-x*100000).round
    new_northing = (n-y*100000).round

    puts "#{grid[y][x]} #{new_easting} #{new_northing}" if @debug

    [grid[y][x], new_easting, new_northing]
  end

  # Iteratively convert from cartesian to lat, long and height.
  def cart_to_llh(x, y, z, ell_type = 'Airy 1830', prec = 10e-9)
    ell      = EllipsoidInfo.new(ell_type)
    a, b, e_2= ell.returnInfo

    # TODO: check quadrant?
    lam      = Math.atan(y.to_f/x.to_f)

    p        = (x**2 + y**2)**0.5
    phi_a    = Math.atan(z.to_f / (p * (1 - e_2)))
    phi_b    = 0

    # Iteratively compute phi.
    while ((phi_a - phi_b).abs > prec)
      v     = a * ((1-e_2 * sin2(phi_a))**-0.5)
      tmp   = phi_a
      phi_a = Math.atan((z + e_2*v*sin(phi_a)) / p)
      phi_b = tmp
    end
    # Compute the height.
    h = p / cos(phi_a) - v

    [phi_a, lam, h]
  end

  # Convert from lat, long and ellipsoid height to cartesian
  def llh_to_cart(phi, lam, h, ell_type = 'Airy 1830')
    ell      = EllipsoidInfo.new(ell_type)
    a, b, e_2= ell.returnInfo

    v        = a * ((1-e_2 * sin2(phi))**-0.5)
    x        = (v+h)*cos(phi) * cos(lam)
    y        = (v+h)*cos(phi) * sin(lam)
    z        = ((1-e_2)*v+h)  * sin(phi)

    [x, y, z]
  end
 
  # ng_to_en operates 5 fig e and n - this function takes a string input, and 
  # converts it to the right sized int. E.g.:
  # 513 becomes 51300, and 00513 stays as 513.
  def fix_ng_input(v)
    if v.length < 5
      m = 10**(5 - v.length)
      v.to_i*m
    else
      v.to_i
    end
  end 
  
  def fix_ng_inputs(e, n)
    [fix_ng_input(e), fix_ng_input(n)]
  end
 
  # Convert from National Grid, back to plain old TM.
  def ng_to_en(prefix, e, n, grid_type = 'GB')
    # E and N multipliers.
    me = 0;
    mn = 0;

    if e.class == String
      e = fix_ng_input(e)
    end
    if n.class == String
      n = fix_ng_input(n)
    end

    g = GridInfo.new(grid_type)
    grid = g.returnInfo

    # Search the grid for the prefix location.
    grid.each_with_index {|row, i|
      row.each_with_index {|entry, j|
        if entry == prefix
          puts "Found #{entry} at #{i} #{j}" if @debug
          mn = i
          me = j
          break
        end
      }
    }

    [e+me*100000, n+mn*100000]  
  end

  def ng_to_ll(prefix, e, n, grid_type = 'GB', proj_type = 'National Grid')
    en, nn =  ng_to_en(prefix, e, n, grid_type)
    # Now convert to lat and long (in radians).
    phi, lam = en_to_ll(en, nn, proj_type)

    [phi, lam]
  end
  
  private :fix_ng_input, :fix_ng_inputs
  
end

# Class for general single coordinate operations.
# Uses the functions from CoordinateConversion
# Implements automatic type conversion. e.g.:
#   a = Coordinate.new(52, 39, 27.2531)
#   a.to_dd(52, 39, 27.2531)
#   => 52.6575703055556
# *Note*: The <tt>to_*</tt> methods do not appear in the docs below because
# they are generated automatically using the listTypes function of 
# CoordinateConversion
class Coordinate
  include CoordinateConversion

  attr_reader :init, :types

  def initialize(*args)
    @dd  = nil
    @dm  = nil
    @dms = nil
    @init = nil

    listTypes

    if args.length == 1
      if args[0].kind_of? Numeric
        @dd = args[0]
        @init = 'dd'
      else 
        raise ArgumentError, "Bad input argument", caller
      end
    elsif args.length == 2
      @dm = args 
      @init = :dm
    elsif args.length == 3
      @dms = args
      @init = 'dms'
    end

    # Build accessor functions.
    @types.each {|type|
      eval "def to_#{type}; get('#{type}'); end"
    }
  end

  def get(type)
    # Return the wanted val, converting if needed.
    if type != @init
      eval "@#{type} = #{@init}_to_#{type}(*@#{@init})"
    else
      eval "@#{@init}"
    end
  end

end

# Parent class for the *Info classes.
class Info
  def initialize(type)
    @type = type
  end
end

# A class containing information about (various) National Grids,
# for use in <tt>ng_to_en</tt> and <tt>en_to_ng</tt>
class GridInfo < Info
  def initialize(type = 'GB')
    super(type)
  end

  # Return info on selected grid.
  def returnInfo
    case @type
    when 'GB'
      grid = [['SV', 'SW', 'SX', 'SY', 'SZ', 'TV', 'TW'],
              ['SQ', 'SR', 'SS', 'ST', 'SU', 'TQ', 'TR'],
              ['SL', 'SM', 'SN', 'SO', 'SP', 'TL', 'TM'], 
              ['SF', 'SG', 'SH', 'SJ', 'SK', 'TF', 'TG'],
              ['SA', 'SB', 'SC', 'SD', 'SE', 'TA', 'TB'],
              ['NV', 'NW', 'NX', 'NY', 'NZ', 'OV', 'OW'],
              ['NQ', 'NR', 'NS', 'NT', 'NU', 'OQ', 'OR'],
              ['NL', 'NM', 'NN', 'NO', 'NP', 'OL', 'OM'],
              ['NF', 'NG', 'NH', 'NJ', 'NK', 'OF', 'OG'],
              ['NA', 'NB', 'NC', 'ND', 'NE', 'OA', 'OB'],
              ['HV', 'HW', 'HX', 'HY', 'HZ', 'JV', 'JW'],
              ['HQ', 'HR', 'HS', 'HT', 'HU', 'JQ', 'JR'],
              ['HL', 'HM', 'HN', 'HO', 'HP', 'JL', 'JM']] 
    else
      raise ArgumentError, 'Bad input argument', caller
    end
    grid
  end
end

# A class containing information about Ellipsiods used in Coordinates.
class EllipsoidInfo < Info
  def initialize(type = 'Airy 1830')
    super(type)
  end

  # Return all of the info about the selected ellipsoid type.
  def returnInfo
    case @type
    when 'Airy 1830' || '1830'
      a = 6377563.396 
      b = 6356256.910
    when 'Airy 1830 modified' || '1890m'
      a = 6377240.189
      b = 6356034.477
    when 'International 1924' || 'Hayford 1909' || 'i1924'
      a = 6378388.000
      b = 6356911.946
    when 'WGS84' || 'GRS80'
      a = 6378137.000
      b = 6356752.3141
    else
      raise ArgumentError, "Bad input argument", caller
    end
    e_2	= (a**2 - b**2)/a**2
    [a, b, e_2]
  end
end

# A class containing information about projections used in Coordinates.
class ProjectionInfo < Info
  def initialize(type = 'National Grid')
    super(type)
  end
  # Return all of the info about the selected projection type.
  def returnInfo
    case @type
    when 'National Grid' || 'OS'
      phi_0 = 49.deg_to_rad 
      lam_0 = -2.deg_to_rad 
      n_0   = -100000
      e_0 	=  400000
      f_0   =  0.9996012717 
      ell = EllipsoidInfo.new('Airy 1830')
    when 'UTM zone 29'
      phi_0 = 0
      lam_0 = -9.deg_to_rad
      n_0   = 0
      e_0   = 500000
      f_0   = 0.9996
      ell = EllipsoidInfo.new('International 1924')
    when 'UTM zone 30'
      phi_0 =  0
      lam_0 = -3.deg_to_rad
      n_0   =  0
      e_0   =  500000
      f_0   =  0.9996
      ell = EllipsoidInfo.new('International 1924')
    when 'UTM zone 31'
      phi_0 =  0
      lam_0 =  3.deg_to_rad
      n_0   =  0
      e_0   =  500000
      f_0   =  0.9996
      ell = EllipsoidInfo.new('International 1924')
    else
      raise ArgumentError, "Bad input argument", caller
    end
    a, b, e_2 = ell.returnInfo
    [phi_0, lam_0, n_0, e_0, f_0, a, b, e_2]
  end
end







