const std = @import("std");

pub fn main() !void {
    const allocator = std.heap.page_allocator;
    const expression = "exp(3,5^2)+3.14*pi*4*(x-6)";
    const tokens = try tokenizeExpression(allocator, expression);
    const rpn = try shuntingYardAlgorithm(allocator, tokens);

    std.debug.print("RPN: ", .{});
    for (rpn) |token| {
        std.debug.print("{s} ", .{token.value});
    }
    std.debug.print("\n", .{});
}

pub fn operandsPerOperator(operator: []const u8) usize {
    if (operator.len == 1)
        switch (operator) {
            '+', '-', '*', '/' => return 2,
            '^' => return 2,
            else => @panic("Invalid RPN"),
        };
}

// pub fn leftBound(f: [][]const u8, operatorIndex: usize) usize {
//     var operandsNeeded = operandsPerOperator(f[operatorIndex][0]);
// }

pub fn derivate(allocator: std.mem.Allocator, f: [][]const u8) [][]const u8 {
    _ = allocator;
    return @constCast(&[_][]const u8{f[1]});
}

test derivate {
    const input: [][]const u8 = @constCast(&[_][]const u8{
        "x",
        "2",
        "*",
    });

    try std.testing.expectEqualSlices([]const u8, &[_][]const u8{"2"}, derivate(std.testing.allocator, input));
}

const TokenType = union(enum) {
    Number: Number,
    Variable: Variable,
    Function: Function,
    Operator: Operator,
    Constants: Constants,
    LeftParen: Paren,
    RightParen: Paren,
};

const Number = struct {
    value: []const u8,
};

const Variable = struct {
    name: []const u8,
};

const Function = struct {
    name: []const u8,
};

const Operator = struct {
    kind: []const u8,
    precedence: u8,
    leftAssociative: bool,
};

const Constants = struct {
    value: []const u8,
};

const Paren = struct {
    left: bool,
};

const Token = struct {
    type: TokenType,
    value: []const u8,
};

const functions = [_]Function{
    .{ .name = "sin" },
    .{ .name = "cos" },
    .{ .name = "tan" },
    .{ .name = "asin" },
    .{ .name = "acos" },
    .{ .name = "atan" },
    .{ .name = "sinh" },
    .{ .name = "cosh" },
    .{ .name = "tanh" },
    .{ .name = "asinh" },
    .{ .name = "acosh" },
    .{ .name = "atanh" },
    .{ .name = "ln" },
    .{ .name = "exp" },
    .{ .name = "sqrt" },
    .{ .name = "abs" },
    .{ .name = "neg" },
    .{ .name = "ceil" },
    .{ .name = "floor" },
    .{ .name = "gamma" },
    .{ .name = "lgamma" },
};

const operators = [_]Operator{
    .{ .kind = "+", .precedence = 1, .leftAssociative = true },
    .{ .kind = "-", .precedence = 1, .leftAssociative = true },
    .{ .kind = "*", .precedence = 2, .leftAssociative = true },
    .{ .kind = "/", .precedence = 2, .leftAssociative = true },
    .{ .kind = "^", .precedence = 3, .leftAssociative = false },
};

const constants = [_]Constants{
    .{ .value = "e" },
    .{ .value = "pi" },
    .{ .value = "phi" },
};

fn isFunction(name: []const u8) bool {
    for (functions) |func| {
        if (std.mem.eql(u8, func.name, name)) {
            return true;
        }
    }
    return false;
}

fn isConstant(name: []const u8) bool {
    for (constants) |constant| {
        if (std.mem.eql(u8, constant.value, name)) {
            return true;
        }
    }
    return false;
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
                defer allocator.free(number);

                var j: usize = 0;
                while (j < end - i) : (j += 1) {
                    number[j] = if (expression[i + j] == ',') '.' else expression[i + j];
                }

                try tokens.append(Token{
                    .type = TokenType{ .Number = .{ .value = try allocator.dupe(u8, number) } },
                    .value = try allocator.dupe(u8, number),
                });
                i = end;
            },
            'a'...'z', 'A'...'Z' => {
                var end = i + 1;
                while (end < expression.len and std.ascii.isAlphabetic(expression[end])) : (end += 1) {}
                const name = expression[i..end];

                const token = if (isFunction(name))
                    Token{ .type = TokenType{ .Function = .{ .name = name } }, .value = name }
                else if (isConstant(name))
                    Token{ .type = TokenType{ .Constants = .{ .value = name } }, .value = name }
                else
                    Token{ .type = TokenType{ .Variable = .{ .name = name } }, .value = name };

                try tokens.append(token);
                i = end;
            },
            '(', ')' => {
                try tokens.append(Token{
                    .type = if (c == '(')
                        TokenType{ .LeftParen = .{ .left = true } }
                    else
                        TokenType{ .RightParen = .{ .left = false } },
                    .value = if (c == '(') "(" else ")",
                });
                i += 1;
            },
            else => {
                const op_str = expression[i .. i + 1];
                for (operators) |op| {
                    if (std.mem.eql(u8, op.kind, op_str)) {
                        try tokens.append(Token{
                            .type = TokenType{ .Operator = op },
                            .value = op_str,
                        });
                        break;
                    }
                }
                i += 1;
            },
        }
    }

    return tokens.toOwnedSlice();
}

pub fn shuntingYardAlgorithm(allocator: std.mem.Allocator, expression: []const Token) ![]const Token {
    if (expression.len == 0) return &[_]Token{};

    var output_queue = try std.ArrayList(Token).initCapacity(allocator, expression.len);
    var operator_stack = try std.ArrayList(Token).initCapacity(allocator, expression.len);
    defer operator_stack.deinit();

    for (expression) |token| {
        switch (token.type) {
            .Number, .Variable, .Constants => try output_queue.append(token),
            .Function, .LeftParen => try operator_stack.append(token),
            .Operator => |current_op| {
                while (operator_stack.items.len > 0) {
                    const top = operator_stack.items[operator_stack.items.len - 1];
                    switch (top.type) {
                        .Operator => |stack_op| {
                            if (stack_op.precedence < current_op.precedence or
                                (stack_op.precedence == current_op.precedence and !current_op.leftAssociative)) break;

                            try output_queue.append(operator_stack.pop());
                        },
                        .LeftParen => break,
                        else => break,
                    }
                }
                try operator_stack.append(token);
            },
            .RightParen => {
                var found_left_paren = false;
                while (operator_stack.items.len > 0) {
                    const top = operator_stack.pop();
                    if (top.type == .LeftParen) {
                        found_left_paren = true;
                        if (operator_stack.items.len > 0 and operator_stack.items[operator_stack.items.len - 1].type == .Function) {
                            try output_queue.append(operator_stack.pop());
                        }
                        break;
                    }
                    try output_queue.append(top);
                }
                if (!found_left_paren) return error.MismatchedParentheses;
            },
        }
    }

    while (operator_stack.items.len > 0) {
        const top = operator_stack.pop();
        if (top.type == .LeftParen or top.type == .RightParen) {
            return error.MismatchedParentheses;
        }
        try output_queue.append(top);
    }

    return output_queue.toOwnedSlice();
}
