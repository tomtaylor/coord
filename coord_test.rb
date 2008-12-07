#! /usr/bin/env ruby

# Tests for coord.rb

require 'coord.rb'
require 'test/unit'

class ConversionTest
  include CoordinateConversion

  attr_reader :types

  def initialize
    listTypes
  end
  
end

module CompTest

  ERR = 1e-6 # Level of error we're willing to accept.
  
  def within_err(v1, v2)
    val = (v1 - v2).abs < ERR
    assert(val, "#{v1} is not close enough to #{v2}")
  end
end

class CoordinatesTest < Test::Unit::TestCase
  include CompTest
  
  def setup
    @en = [651409.902938192, 313177.270343057]
    @en_i = [651410, 313177]
    @ng = ["TG", 51410, 13177]
    @ng_s = ["TG", '51410', '13177']
    @cart = [3874938.84930327, 116218.623758513, 5047168.20839729]
    @h = 24.7000000048429
    coord = Coordinate.new
    @phi = coord.dms_to_dd(52, 39, 27.2531).deg_to_rad
    @lam = coord.dms_to_dd(1, 43, 4.5177).deg_to_rad
  end
  
  def test_coordinates_ll_to_en
    test =  Coordinates.new
    o = test.ll_to_en(@phi, @lam)
    
    o.each_with_index {|v, i|
      within_err(v, @en[i])
    }
  end
 
  def test_coordinates_en_to_ll
    test =  Coordinates.new
    o = test.en_to_ll(*@en)
    o.each_with_index {|v, i|
      within_err(v, [@phi, @lam][i])
    }
  end

  def test_coordinates_llh_to_cart
    test = Coordinates.new
    o = test.llh_to_cart(@phi, @lam, @h)

    o.each_with_index {|v, i|
      within_err(v, @cart[i])
    }
  end
 
  def test_coordinates_cart_to_llh
    test = Coordinates.new
    o = test.cart_to_llh(*@cart)
    # TODO: Replace llh_rad with [phi, lam, h]
    
    o.each_with_index {|v, i|
      within_err(v, [@phi, @lam, @h][i])
    }
  end
 
  def test_coordinates_en_to_ng
    test = Coordinates.new
    val = test.en_to_ng(*@en)
    assert(val == @ng)
  end
  
  def test_coordinates_ng_to_en
    test = Coordinates.new
    val = test.ng_to_en(*@ng) 
    # Test against integer E and N.
    assert(val == @en_i)

    # Now test string handling:
    val = test.ng_to_en(*@ng_s) 
    # Test against integer E and N.
    assert(val == @en_i)

    #..and now test 6 fig. handling.
    val = test.ng_to_en(@ng_s[0], @ng_s[1][0..2], @ng_s[2][0..2]) 
    # Test against integer E and N.
    en_t = @en_i.collect {|v| (((v/100).floor)*100)}
    assert(val == en_t)
  end
end


class CoordinateTest < Test::Unit::TestCase
  include CompTest
  
  def setup
    @types = {'dms' => [45, 22, 38], 
              'dd'  => 45.3772222222222, 
              'dm'  => [45, 22.63333333]}
  end

  def test_coordinate_types
    # Test to see if all inputs can be handled.
      @types.each {|k, v|
      test = Coordinate.new(*v)
      assert(k, test.init)
    }

    begin   
      test_bad = Coordinate.new('foo')
    rescue
      assert($!.class, ArgumentError)
    end
    
  end

  def test_conversion_module
    # Test the conversions between various coord types.
    # Ignore rounding error.
    test = ConversionTest.new
    val = false

    @types.each {|k, v|
      assert(test.types.include?(k))
      @types.each{|k2, v2|
        if v.class == Array
          v = v.join(', ')
        end
        ret = eval "test.#{k}_to_#{k2}(*[#{v}])"
        if ret.class == Array
          ret.each_with_index {|s, i|
            within_err(s, v2[i])
          }
        else
          within_err(ret, v2)
        end
      }
    }
  end
end
