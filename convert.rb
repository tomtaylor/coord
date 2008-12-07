#! /usr/bin/env ruby

# :Title: convert.rb: Convert sets of coordinates on the command line
# Author: Matt Foster (mailto:mattfoster@clara.co.uk)
# Copyright: Matt Foster (mailto:mattfoster@clara.co.uk)
# Convert coordinate types on the command line.
# For help, try convert.rb -h
# E.g.:
#   ruby convert.rb -i 'TG 51410 13177' -o dms 
#   E 1D 43M 4.52210295347086S
#   N 52D 39M 27.2442040793669S


require 'coord.rb'
require 'optparse' 
require 'ostruct'

def split_coord(c)
  if c.length % 2 != 0
    puts 'Error expected an even length coord set.'
    puts p
    exit(1)
  elsif c.class != String
    puts 'Expected string input.'
    puts p
    exit(1)
  end
   
  # Split the input coordinate in two.
  e = c[0..c.length/2-1]
  n = c[c.length/2..-1]

  [e.to_i, n.to_i]
end

def get_ll_prefix(phi, lam)
  if lam < 0
    lam_prefix = 'W'
  else
    lam_prefix = 'E'
  end
  if phi < 0
    phi_prefix = 'W'
  else
    phi_prefix = 'N'
  end
  [phi_prefix, lam_prefix]
end

def put_dms(phi, lam)
  # Now prettify the output.
  if phi.class != Coordinate and lam.class != Coordinate
    phi_o = Coordinate.new(phi.rad_to_deg).to_dms
    lam_o = Coordinate.new(lam.rad_to_deg).to_dms
  else
    phi_o = phi.to_dms
    lam_o = lam.to_dms
  end
  phi_prefix, lam_prefix = get_ll_prefix(phi_o[0], lam_o[0])

  puts "#{lam_prefix} #{lam_o[0].abs}D #{lam_o[1]}M #{lam_o[2]}S"
  puts "#{phi_prefix} #{phi_o[0].abs}D #{phi_o[1]}M #{phi_o[2]}S"
end

def put_dm(phi, lam)
  # Now prettify the output.
  lam_o = Coordinate.new(lam.rad_to_deg).to_dm
  phi_o = Coordinate.new(phi.rad_to_deg).to_dm

  phi_prefix, lam_prefix = get_ll_prefix(phi, lam)

  puts "#{phi_prefix} #{phi_o[0].abs}D #{phi_o[1]}M"
  puts "#{lam_prefix} #{lam_o[0].abs}D #{lam_o[1]}M"
end

def put_dd(phi, lam)
  # Now prettify the output.
  lam_o = Coordinate.new(lam.rad_to_deg).to_dd
  phi_o = Coordinate.new(phi.rad_to_deg).to_dd

  phi_prefix, lam_prefix = get_ll_prefix(phi, lam)

  puts "#{phi_prefix} #{phi_o.abs}D"
  puts "#{lam_prefix} #{lam_o.abs}D"
end

def put_ng(p, ee, nn)
  puts "#{p} #{ee} #{nn}"
end

def put_en(e, n)
  puts "E #{e}, N #{n}"
end

# Parse the command line arguments. Users OptionParser.
opt = OpenStruct.new
# Default values.
opt.verbose = false
# Nice simple way of parsing opts.
p = OptionParser.new {|p|
  p.on('-i input_coords', '--input', 'Input Coordinates', String) {|o| 
    opt.input = o unless o.empty?
  }
  p.on('-o [TYPE]', '--output', [:dms, :dd, :dm, :en, :ng], 
    "Select output type (dms, dd, dm, en, ng)")   {|o|
    opt.output = o
  }
  p.on('-v', '--[no-]verbose', 'Run verbosely') {|o|
    opt.verbose = o
  }
  p.on_tail("-h", "--help", "Show this message") {
    puts p
    puts "Convert coordinates on the command line."
    exit(1)
  }
}
p.parse!

# If there aren't any arguments, show help.
if opt.input == nil or opt.output == nil
  puts p
  exit(1)
end

c  = Coordinates.new

# get args for necessary conversion.

case opt.input
when /^(\d+)\s+(\d+)$/
  #puts 'EN'
  input = :en
  e = $1.to_i
  n = $2.to_i
  args = [e, n]

when /^(\d+)$/
  #puts 'EN: No spaces'
  input = :en
  e, n = split_coord($1)
  args = [e, n]

when /^(\w\w)\s*(\d+)\s+(\d+)$/
  #puts 'National Grid Format'
  input = :ng
  args = [$1, $2.to_i, $3.to_i]

when /^(\w\w)\s*(\d+)$/
  #puts 'NG: no spaces.'
  input = :ng
  p = $1
  e, n = split_coord($2)
  args = [p, e, n]

when /^(\w)\s*(\d+)\w? (\d+\.?\d+)\w? (\d+\.?\d+)\w?, (\w)\s*(\d+)\w? (\d+\.?\d+)\w? (\d+\.?\d+)\w?$/
  #puts 'DMS'
  input = :ll
  # Todo resolve leading NSEW to +/-
  de = $2.to_i
  if $1.upcase == 'S'
    de = de * -1
  end
 
  dn = $6.to_i
  if $5.upcase == 'W'
    dn = dn * -1
  end
  
  c1 = Coordinate.new(de, $3.to_f, $4.to_f)
  c2 = Coordinate.new(dn, $7.to_f, $8.to_f)

  args = [c1, c2]
 
when /^(\w)\s*(\d+)\w (\d+\.\d+)\w$/
  puts 'DM'
  puts $1, $2, $3, $4
 
when /^(\w)\s*(\d+\.\d+)\w$/
  puts 'D'
  puts $1, $2, $3, $4
 
end

if opt.output == :dms or opt.output == :dd or opt.output == :dm
  op = :ll
else
  op = opt.output
end

if input == op
  conv = "put_#{opt.output}(*args)"
else
  o = eval("c.#{input}_to_#{op}(*args)")
  conv =  "put_#{opt.output}(*o)"
end

eval(conv)
