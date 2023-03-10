function tests = sumsumtest
tests = functiontests(localfunctions);
end

function testsumsum(testcase)

example = [1 2];
expectation = 2
verifyEqual(testcase,example,expectation)

end
