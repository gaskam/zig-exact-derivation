const std = @import("std");

pub fn main() !void {
    const allocator = std.heap.page_allocator;
    const expression = "exp(3,5^2 * x)+3.14*pi*4*(x-6)";
    const tokens = try tokenizeExpression(allocator, expression);
    const rpn = try shuntingYardAlgorithm(allocator, tokens);

    std.debug.print("RPN: ", .{});
    for (rpn) |token| {
        std.debug.print("{s} ", .{try tagToStr(allocator, token)});
    }
    std.debug.print("\n", .{});
}

/// Returns the derivative of the given expression
/// as an array of tokens in reverse polish notation.
pub fn getDerivativeTokens(expression: []const u8) ![]Token {
    const allocator = std.heap.page_allocator;
    const tokens = try tokenizeExpression(allocator, expression);
    const rpn = try shuntingYardAlgorithm(allocator, tokens);

    var output = std.ArrayList(Token).init(allocator);
    try derivate(&output, rpn);

    return output.toOwnedSlice();
}

/// Calculates the derivative of the given expression at the given x value.
pub fn calcDerivative(expression: []const u8, x: f64) !f64 {
    const allocator = std.heap.page_allocator;
    const tokens = try getDerivativeTokens(expression);
    const result = try evaluateExpression(allocator, tokens, x);
    return result;
}

/// Calculates the derivative of the given expression at the given x value
/// with the given precision. Maximum precision is 1e-50.
pub fn calcDerivativePrecise(expression: []const u8, x: f64, precision: f64) !f64 {
    _ = .{expression, x, precision};
    //TODO: Implement    
}

pub fn tagToStr(allocator: std.mem.Allocator, token: Token) ![]const u8 {
    switch (token) {
        .Number => |number| {
            switch (number) {
                .integer => return try std.fmt.allocPrint(allocator, "{d}", .{number.integer}),
                .decimal => return number.decimal,
            }
        },
        .Variable => |variable| return variable.name,
        .Function => |function| return @tagName(function),
        .Operator => |operator| return @tagName(operator),
        .Constants => |constants| return @tagName(constants),
        .Parenthesis => return "Parenthesis",
        .Null => return "Null",
    }
}

const Token = union(enum) {
    Number: Number,
    // This should only be reached in/after the derivation step
    Null: void,
    Constants: Constant,
    Variable: Variable,
    Function: Function,
    Operator: Operator,
    Parenthesis: Parenthesis,
};

const Number = union(enum) {
    integer: i64,
    decimal: []const u8,
};

const Variable = struct {
    name: []const u8,
};

const Parenthesis = enum {
    left,
    right,
};

const Function = enum {
    sin,
    cos,
    tan,
    asin,
    acos,
    atan,
    sinh,
    cosh,
    tanh,
    asinh,
    acosh,
    atanh,
    ln,
    exp,
    sqrt,
    neg,
};

const Operator = enum {
    @"+",
    @"-",
    @"*",
    @"/",
    @"^",
};

const OperatorInfo = struct {
    precedence: u8,
    leftAssociative: bool,
};

const operators = struct {
    const @"+": OperatorInfo = .{ .precedence = 1, .leftAssociative = true };
    const @"-": OperatorInfo = .{ .precedence = 1, .leftAssociative = true };
    const @"*": OperatorInfo = .{ .precedence = 2, .leftAssociative = true };
    const @"/": OperatorInfo = .{ .precedence = 2, .leftAssociative = true };
    const @"^": OperatorInfo = .{ .precedence = 3, .leftAssociative = false };

    fn getOperatorInfo(op: Operator) OperatorInfo {
        return switch (op) {
            .@"+" => operators.@"+",
            .@"-" => operators.@"-",
            .@"*" => operators.@"*",
            .@"/" => operators.@"/",
            .@"^" => operators.@"^",
        };
    }
};

const Constant = enum {
    e,
    pi,
    phi,
};

pub fn getFunction(name: []const u8) ?Function {
    return std.meta.stringToEnum(Function, name);
}

pub fn getConstant(name: []const u8) ?Constant {
    return std.meta.stringToEnum(Constant, name);
}

pub fn getOperator(name: []const u8) ?Operator {
    return std.meta.stringToEnum(Operator, name);
}

fn tokenizeExpression(allocator: std.mem.Allocator, expression: []const u8) ![]Token {
    var tokens = std.ArrayList(Token).init(allocator);
    var i: usize = 0;

    while (i < expression.len) {
        const c = expression[i];

        switch (c) {
            ' ', '\t', '\n', '\r' => i += 1,
            '0'...'9' => {
                var end = i + 1;
                var has_decimal = false;
                while (end < expression.len) : (end += 1) {
                    const next = expression[end];
                    if (!std.ascii.isDigit(next) and next != '.' and next != ',') break;
                    if (next == '.' or next == ',') {
                        if (has_decimal) return error.InvalidNumber;
                        has_decimal = true;
                    }
                }

                var number = try allocator.alloc(u8, end - i);

                var j: usize = 0;
                while (j < end - i) : (j += 1) {
                    number[j] = if (expression[i + j] == ',') '.' else expression[i + j];
                }

                if (has_decimal) {
                    try tokens.append(Token{ .Number = .{ .decimal = number } });
                } else {
                    try tokens.append(Token{ .Number = .{ .integer = try std.fmt.parseInt(i64, number, 10) } });
                    allocator.free(number);
                }

                i = end;
            },
            'a'...'z', 'A'...'Z' => {
                var end = i + 1;
                while (end < expression.len and std.ascii.isAlphabetic(expression[end])) : (end += 1) {}
                const name = expression[i..end];

                const token = if (getFunction(name)) |function|
                    Token{ .Function = function }
                else if (getConstant(name)) |constant|
                    Token{ .Constants = constant }
                else
                    Token{ .Variable = .{ .name = name } };

                try tokens.append(token);
                i = end;
            },
            '(', ')' => {
                try tokens.append(Token{
                    .Parenthesis = if (c == '(') Parenthesis.left else Parenthesis.right,
                });
                i += 1;
            },
            '+', '-', '*', '/', '^' => {
                try tokens.append(Token{
                    .Operator = getOperator(expression[i .. i + 1]) orelse unreachable,
                });
                i += 1;
            },
            else => return error.InvalidCharacter,
        }
    }

    return tokens.toOwnedSlice();
}

pub fn shuntingYardAlgorithm(allocator: std.mem.Allocator, expression: []const Token) ![]Token {
    if (expression.len == 0) return &[_]Token{};

    var output_queue = try std.ArrayList(Token).initCapacity(allocator, expression.len);
    var operator_stack = try std.ArrayList(Token).initCapacity(allocator, expression.len);
    defer operator_stack.deinit();

    for (expression) |token| {
        switch (token) {
            .Number, .Variable, .Constants => try output_queue.append(token),
            .Function => try operator_stack.append(token),
            .Operator => |current_op| {
                while (operator_stack.items.len > 0) {
                    const top = operator_stack.items[operator_stack.items.len - 1];
                    switch (top) {
                        .Operator => |stack_op| {
                            const stack_op_info = operators.getOperatorInfo(stack_op);
                            const current_op_info = operators.getOperatorInfo(current_op);
                            if (stack_op_info.precedence < current_op_info.precedence or
                                (stack_op_info.precedence == current_op_info.precedence and !current_op_info.leftAssociative)) break;

                            try output_queue.append(operator_stack.pop());
                        },
                        else => break,
                    }
                }
                try operator_stack.append(token);
            },
            .Parenthesis => |parenthesis| {
                if (parenthesis == Parenthesis.left) {
                    try operator_stack.append(token);
                    continue;
                }
                var found_left_paren = false;
                while (operator_stack.items.len > 0) {
                    const top = operator_stack.pop();
                    if (top == .Parenthesis) {
                        found_left_paren = true;
                        if (operator_stack.items.len > 0 and operator_stack.items[operator_stack.items.len - 1] == .Function) {
                            try output_queue.append(operator_stack.pop());
                        }
                        break;
                    }
                    try output_queue.append(top);
                }
                if (!found_left_paren) return error.MismatchedParentheses;
            },
            .Null => unreachable,
        }
    }

    while (operator_stack.items.len > 0) {
        const top = operator_stack.pop();
        if (top == .Parenthesis) {
            return error.MismatchedParentheses;
        }
        try output_queue.append(top);
    }

    return output_queue.toOwnedSlice();
}

pub fn evaluateExpression(allocator: std.mem.Allocator, tokens: []const Token, x: f64) !f64 {
    var stack = std.ArrayList(f64).init(allocator);
    defer stack.deinit();

    for (tokens) |token| {
        switch (token) {
            .Number => |number| {
                const value = switch (number) {
                    .integer => |i| @as(f64, @floatFromInt(i)),
                    .decimal => |d| try std.fmt.parseFloat(f64, d),
                };
                try stack.append(value);
            },
            .Variable => try stack.append(x),
            .Constants => |constant| {
                const value: f64 = switch (constant) {
                    .e => std.math.e,
                    .pi => std.math.pi,
                    .phi => std.math.phi,
                };
                try stack.append(value);
            },
            .Function => |function| {
                if (stack.items.len < 1) return error.InvalidExpression;
                const arg = stack.pop();
                const result = switch (function) {
                    .sin => @sin(arg),
                    .cos => @cos(arg),
                    .tan => @tan(arg),
                    .asin => std.math.asin(arg),
                    .acos => std.math.acos(arg),
                    .atan => std.math.atan(arg),
                    .sinh => std.math.sinh(arg),
                    .cosh => std.math.cosh(arg),
                    .tanh => std.math.tanh(arg),
                    .asinh => std.math.asinh(arg),
                    .acosh => std.math.acosh(arg),
                    .atanh => std.math.atanh(arg),
                    .ln => @log(arg),
                    .exp => @exp(arg),
                    .sqrt => @sqrt(arg),
                    .neg => -arg,
                };
                try stack.append(result);
            },
            .Operator => |operator| {
                if (stack.items.len < 2) return error.InvalidExpression;
                const b = stack.pop();
                const a = stack.pop();
                const result = switch (operator) {
                    .@"+" => a + b,
                    .@"-" => a - b,
                    .@"*" => a * b,
                    .@"/" => a / b,
                    .@"^" => std.math.pow(f64, a, b),
                };
                try stack.append(result);
            },
            .Parenthesis => return error.UnexpectedToken,
            .Null => return error.UnexpectedToken,
        }
    }

    if (stack.items.len != 1) return error.InvalidExpression;
    return stack.items[0];
}

const Range = struct {
    start: usize,
    end: usize,
};

pub fn getLeftMostOperand(f: []Token) usize {
    var needed: usize = 1;
    var i: usize = f.len - 1;
    while (needed > 0) : ({
        needed -= 1;
        i -= 1;
    }) {
        switch (f[i]) {
            .Number, .Null, .Variable, .Constants => {},
            .Function => needed += 1,
            .Operator => needed += 2,
            .Parenthesis => unreachable,
        }
    }
    return i + 1;
}

pub fn derivate(output: *std.ArrayList(Token), f: []Token) !void {
    switch (f[f.len - 1]) {
        .Number => try output.append(Token{ .Null = {} }),
        .Null => unreachable,
        .Constants => try output.append(Token{ .Null = {} }),
        .Variable => try output.append(Token{ .Number = .{ .integer = 1 } }),
        .Function => |function| {
            switch (function) {
                .neg => {
                    try derivate(output, f[0 .. f.len - 1]);
                    try output.append(Token{ .Function = .neg });
                },
                else => @panic("Not implemented"),
            }
        },
        .Operator => |operator| {
            const pivot = getLeftMostOperand(f[0 .. f.len - 1]);
            const u = f[0..pivot];
            const v = f[pivot .. f.len - 1];
            var changed: u2 = 0;
            switch (operator) {
                .@"+" => {
                    try derivate(output, u);
                    if (output.items[output.items.len - 1] == .Null) {
                        _ = output.pop();
                    } else changed += 1;

                    try derivate(output, v);
                    if (output.items[output.items.len - 1] == .Null) {
                        _ = output.pop();
                    } else changed += 1;

                    if (changed == 2) try output.append(Token{ .Operator = .@"+" });
                    if (changed == 0) try output.append(Token{ .Null = {} });
                },
                .@"-" => {
                    try derivate(output, u);
                    if (output.items[output.items.len - 1] == .Null) {
                        _ = output.pop();
                    } else changed += 1;

                    try derivate(output, v);
                    if (output.items[output.items.len - 1] == .Null) {
                        _ = output.pop();
                    } else {
                        if (changed == 0) try output.append(Token{ .Function = .neg });
                        changed += 1;
                    }

                    if (changed == 2) try output.append(Token{ .Operator = .@"-" });
                    if (changed == 0) try output.append(Token{ .Null = {} });
                },
                .@"*" => {
                    try derivate(output, u);
                    if (output.items[output.items.len - 1] == .Null) {
                        _ = output.pop();
                    } else {
                        try output.appendSlice(v);
                        try output.append(Token{ .Operator = .@"*" });
                        changed += 1;
                    }

                    try derivate(output, v);
                    if (output.items[output.items.len - 1] == .Null) {
                        _ = output.pop();
                    } else {
                        try output.appendSlice(u);
                        try output.append(Token{ .Operator = .@"*" });
                        changed += 1;
                    }

                    if (changed == 2) try output.append(Token{ .Operator = .@"+" });
                    if (changed == 0) try output.append(Token{ .Null = {} });
                },
                .@"/" => {
                    try derivate(output, u);
                    if (output.items[output.items.len - 1] == .Null) {
                        _ = output.pop();
                    } else {
                        try output.appendSlice(v);
                        try output.append(Token{ .Operator = .@"*" });
                        changed += 1;
                    }

                    try derivate(output, v);
                    if (output.items[output.items.len - 1] == .Null) {
                        _ = output.pop();
                    } else {
                        try output.appendSlice(u);
                        try output.append(Token{ .Operator = .@"*" });
                        if (changed == 0) try output.append(Token{ .Function = .neg });
                        changed += 1;
                    }

                    if (changed == 2) try output.append(Token{ .Operator = .@"-" });
                    if (changed == 0) try output.append(Token{ .Null = {} });

                    try output.appendSlice(v);
                    try output.appendSlice(&[_]Token{ Token{ .Operator = .@"^" }, Token{ .Number = .{ .integer = 2 } }, Token{ .Operator = .@"/" } });
                },
                .@"^" => {
                    // u * v
                    try output.appendSlice(u);
                    try output.appendSlice(v);
                    try output.append(Token{ .Operator = .@"*" });

                    // v' ln (u)
                    try derivate(output, v);
                    if (output.items[output.items.len - 1] == .Null) {
                        _ = output.pop();
                    } else {
                        try output.appendSlice(u);
                        try output.append(Token{ .Function = .ln });
                        try output.append(Token{ .Operator = .@"*" });
                        changed += 1;
                    }

                    // u' * v / u
                    try derivate(output, u);
                    if (output.items[output.items.len - 1] == .Null) {
                        _ = output.pop();
                    } else {
                        try output.appendSlice(v);
                        try output.append(Token{ .Operator = .@"*" });
                        try output.appendSlice(u);
                        try output.append(Token{ .Operator = .@"/" });
                        changed += 1;
                    }
                },
            }
        },
        else => {},
    }
}

test derivate {
    const allocator = std.heap.page_allocator;
    const expression = "3,5+2*x+3.14*pi*4*(x-6)";
    const output = try getDerivativeTokens(expression);

    std.debug.print("RPN Derivative: ", .{});
    for (output) |token| {
        std.debug.print("{s} ", .{try tagToStr(allocator, token)});
    }
    std.debug.print("\n", .{});

    const x = 3.0;
    const result = try calcDerivative(expression, x);
    std.debug.print("Result: {d}\n", .{result});
}
